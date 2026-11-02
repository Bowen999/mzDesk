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
up=0
for i in $(seq 1 20); do
  if curl -fsS -o /dev/null http://127.0.0.1:8599/api/health; then echo "   mzdesk is up ($(systemctl is-active mzdesk))"; up=1; break; fi
  sleep 1
done
if [ "$up" = 0 ]; then echo "!! mzdesk did not answer in 20 s"; journalctl -u mzdesk -n 30 --no-pager; exit 1; fi

# Release for the local version (staged by deploy.sh). Only now, so version.json never points at a site that is not up.
# root owns the files and the app (user mzdesk) only reads them: a taken-over app cannot swap the wheel people install.
# Wheel and checksum first, version.json last and in one step (mv), so a reader sees the old release or the new one, never half of it.
REL=/var/lib/mzdesk/releases
if [ -f "$STAGE/release/version.json" ]; then
  echo "   release $(sed -n 's/.*"version": *"\([^"]*\)".*/\1/p' "$STAGE/release/version.json") -> $REL"
  install -d -o root -g mzdesk -m 750 "$REL"
  for f in "$STAGE"/release/mzdesk-*-py3-none-any.whl "$STAGE"/release/mzdesk-*-py3-none-any.whl.sha256; do
    install -o root -g mzdesk -m 640 "$f" "$REL/"
  done
  install -o root -g mzdesk -m 640 "$STAGE/release/version.json" "$REL/.version.json.new"
  mv -f "$REL/.version.json.new" "$REL/version.json"
  # keep the newest five wheels (and their checksum files)
  ls -1 "$REL"/mzdesk-*-py3-none-any.whl | sort -V | head -n -5 | while read -r old; do rm -f "$old" "$old.sha256"; done
  rm -rf "$STAGE/release"
  curl -fsS -o /dev/null http://127.0.0.1:8599/version.json && echo "   /version.json is served"
fi
