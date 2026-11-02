#!/bin/bash
# One-time move of the public address to https://mzdesk.com/ : certificate, the new build, the nginx site, and a redirect from the
# old https://bowenyang.org/mzdesk/ . Run it after the DNS records point at the server (A mzdesk.com and www -> 143.198.233.45).
# Safe to run again. Every changed nginx file is backed up to /root/nginx-backup on the server, and nginx -t guards each step.
set -euo pipefail
cd "$(dirname "$0")/.."
IP=143.198.233.45
HOST=root@$IP
SSH="ssh -i $HOME/.ssh/mzdesk_deploy -o IdentitiesOnly=yes"

echo "== DNS"
for name in mzdesk.com www.mzdesk.com; do
  got=$(dig +short "$name" A | grep -E '^[0-9.]+$' | sort | tr '\n' ' ')
  echo "   $name -> ${got:-nothing}"
  if [ "${got// /}" != "$IP" ]; then echo "!! $name must resolve to $IP only (Porkbun: DNS records). Stopping."; exit 1; fi
done

echo "== certificate"
$SSH "$HOST" 'mkdir -p /root/mzdesk-deploy /var/www/certbot'
rsync -a -e "$SSH" deploy/nginx-mzdesk-http.conf deploy/00-mzdesk-limits.conf "$HOST:/root/mzdesk-deploy/"
$SSH "$HOST" 'bash -s' <<'REMOTE'
set -euo pipefail
S=/root/mzdesk-deploy; B=/root/nginx-backup; TS=$(date +%Y%m%d-%H%M%S); mkdir -p $B
if [ -f /etc/letsencrypt/live/mzdesk.com/fullchain.pem ]; then echo "   already there"; exit 0; fi
install -m 644 $S/00-mzdesk-limits.conf /etc/nginx/conf.d/00-mzdesk-limits.conf
[ -f /etc/nginx/conf.d/mzdesk.conf ] && cp -p /etc/nginx/conf.d/mzdesk.conf $B/mzdesk.conf.$TS
install -m 644 $S/nginx-mzdesk-http.conf /etc/nginx/conf.d/mzdesk.conf
nginx -t && systemctl reload nginx
if ls -d /etc/letsencrypt/accounts/acme-v02.api.letsencrypt.org/directory/* >/dev/null 2>&1; then MAIL=""; else MAIL="-m feedback@example.org"; fi
certbot certonly --webroot -w /var/www/certbot -d mzdesk.com -d www.mzdesk.com --non-interactive --agree-tos $MAIL
REMOTE

echo "== deploy (code, unit, https site)"
deploy/deploy.sh

echo "== bowenyang.org/mzdesk/ -> https://mzdesk.com/"
$SSH "$HOST" 'bash -s' <<'REMOTE'
set -euo pipefail
B=/root/nginx-backup; TS=$(date +%Y%m%d-%H%M%S); mkdir -p $B
cp -p /etc/nginx/conf.d/bowenyang.conf $B/bowenyang.conf.$TS
python3 - <<'PY'
p = "/etc/nginx/conf.d/bowenyang.conf"
s = open(p).read()
if "https://mzdesk.com" in s:
    print("   already redirects")
else:
    start = s.index("    # mzDesk lives under /mzdesk/")
    end = s.index("    # Old root-level compound pages")
    s = s[:start] + """    # mzDesk moved to https://mzdesk.com/ ; old links keep working (the path after /mzdesk/ is kept)
    location = /mzdesk { return 301 https://mzdesk.com/; }
    location /mzdesk/ { rewrite ^/mzdesk/(.*)$ https://mzdesk.com/$1 permanent; }

""" + s[end:]
    open(p, "w").write(s)
    print("   bowenyang.conf: /mzdesk/ now redirects")
PY
if nginx -t; then systemctl reload nginx; else
  echo "!! nginx -t failed: restoring bowenyang.conf"; cp -p $B/bowenyang.conf.$TS /etc/nginx/conf.d/bowenyang.conf; nginx -t && systemctl reload nginx; exit 1
fi
REMOTE

echo "== check"
curl -s -o /dev/null -w "   https://mzdesk.com/api/health           -> %{http_code}\n" https://mzdesk.com/api/health
curl -s -o /dev/null -w "   https://www.mzdesk.com/                  -> %{http_code} %{redirect_url}\n" https://www.mzdesk.com/
curl -s -o /dev/null -w "   http://mzdesk.com/                       -> %{http_code} %{redirect_url}\n" http://mzdesk.com/
curl -s -o /dev/null -w "   https://bowenyang.org/mzdesk/modules/x   -> %{http_code} %{redirect_url}\n" https://bowenyang.org/mzdesk/modules/x
