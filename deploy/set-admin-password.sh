#!/bin/bash
# Set the admin password of https://mzdesk.com/stats and /reports (one password for both). You type it here, hidden and twice; it is turned
# into a salted hash on this Mac, only the hash goes to /etc/mzdesk/mzdesk.env on the server, and the service restarts. Then the script
# checks both pages with the new password and with a wrong one. The password itself is never written to a file or a command line.
set -euo pipefail
cd "$(dirname "$0")/.."
HOST=root@143.198.233.45
SSH="ssh -i $HOME/.ssh/mzdesk_deploy -o IdentitiesOnly=yes"
PY=${PYTHON:-python3}

read -rsp "New admin password: " pw; echo
read -rsp "Again: " again; echo
[ "$pw" = "$again" ] || { echo "!! The two passwords differ."; exit 1; }
lines=$(printf '%s' "$pw" | $PY -m mzdesk.api.admin --stdin)

echo "== server"
printf '%s\n' "$lines" | $SSH "$HOST" 'f=/etc/mzdesk/mzdesk.env; t=$(mktemp); { grep -v "^MZDESK_ADMIN_\(SALT\|HASH\)=" "$f" || true; cat; } > "$t"; install -m 600 "$t" "$f"; rm -f "$t"; systemctl restart mzdesk; for i in $(seq 1 20); do curl -fsS -o /dev/null http://127.0.0.1:8599/api/health && exit 0; sleep 1; done; echo "!! mzdesk did not come back"; exit 1'

echo "== check"
body() { printf '%s' "$1" | $PY -c 'import sys, json; print(json.dumps({"password": sys.stdin.read()}))'; }
for page in reports stats; do
  ok=$(body "$pw" | curl -s -o /dev/null -w "%{http_code}" -X POST "https://mzdesk.com/api/$page" -H 'content-type: application/json' -d @-)
  bad=$(body "not-$pw" | curl -s -o /dev/null -w "%{http_code}" -X POST "https://mzdesk.com/api/$page" -H 'content-type: application/json' -d @-)
  echo "   /api/$page: new password -> $ok (want 200), wrong password -> $bad (want 401)"
  [ "$ok" = 200 ] && [ "$bad" = 401 ] || { echo "!! the check failed"; exit 1; }
done
if [ -f "$HOME/.mzdesk/admin-password.txt" ]; then rm -f "$HOME/.mzdesk/admin-password.txt"; echo "   removed the old generated password file (~/.mzdesk/admin-password.txt)"; fi
echo "== done: /stats and /reports now use the password you typed"
