#!/bin/bash
# Runs on the server (deploy.sh pipes it over ssh): installs the files staged in /root/mzdesk-deploy, tests nginx (and puts the old
# files back if the test fails), restarts the service and checks that it answers. Safe to run again.
set -euo pipefail
STAGE=/root/mzdesk-deploy
BACKUP=/root/nginx-backup
TS=$(date +%Y%m%d-%H%M%S)
mkdir -p "$BACKUP"

# install <staged file> <target>: copy when different; the old one is kept in $BACKUP. Prints "changed" when it copied.
install_if_changed() {
  local src=$1 dst=$2
  if [ -f "$dst" ] && cmp -s "$src" "$dst"; then return 1; fi
  [ -f "$dst" ] && cp -p "$dst" "$BACKUP/$(basename "$dst").$TS"
  install -m 644 "$src" "$dst"
  echo "   installed $dst"
}

unit_changed=0; nginx_changed=0
install_if_changed "$STAGE/mzdesk.service" /etc/systemd/system/mzdesk.service && unit_changed=1 || true
install_if_changed "$STAGE/00-mzdesk-limits.conf" /etc/nginx/conf.d/00-mzdesk-limits.conf && nginx_changed=1 || true
if [ -f /etc/letsencrypt/live/mzdesk.com/fullchain.pem ]; then
  install_if_changed "$STAGE/nginx-mzdesk.conf" /etc/nginx/conf.d/mzdesk.conf && nginx_changed=1 || true
else
  echo "   no certificate for mzdesk.com yet: nginx site not installed (run deploy/cutover.sh once DNS points here)"
fi

if [ "$nginx_changed" = 1 ]; then
  if nginx -t 2>/dev/null; then
    systemctl reload nginx && echo "   nginx reloaded"
  else
    nginx -t || true
    echo "!! nginx -t failed; putting the previous files back"
    for f in mzdesk.service 00-mzdesk-limits.conf mzdesk.conf; do
      old=$(ls -1t "$BACKUP/$f".* 2>/dev/null | head -1 || true)
      case $f in mzdesk.service) dst=/etc/systemd/system/$f;; *) dst=/etc/nginx/conf.d/$f;; esac
      [ -n "$old" ] && [ "$(stat -c %Y "$old")" -ge "$(date -d '-5 minutes' +%s)" ] && cp -p "$old" "$dst" && echo "   restored $dst"
    done
    nginx -t && systemctl reload nginx
    exit 1
  fi
fi

[ "$unit_changed" = 1 ] && systemctl daemon-reload
echo "   python packages"
/opt/mzdesk/.venv/bin/python -m pip install -q --disable-pip-version-check -r "$STAGE/requirements-server.txt"
chown -R root:mzdesk /opt/mzdesk && chmod -R u=rwX,g=rX,o= /opt/mzdesk
systemctl restart mzdesk
for i in $(seq 1 20); do
  if curl -fsS -o /dev/null http://127.0.0.1:8599/api/health; then echo "   mzdesk is up ($(systemctl is-active mzdesk))"; exit 0; fi
  sleep 1
done
echo "!! mzdesk did not answer in 20 s"; journalctl -u mzdesk -n 30 --no-pager; exit 1
