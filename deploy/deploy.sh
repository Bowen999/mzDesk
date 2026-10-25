#!/bin/bash
# Update the live site (https://bowenyang.org/mzdesk/) from this checkout: build the frontend for /mzdesk/, sync code, restart.
# Models / data under /var/lib/mzdesk are not touched (first-time setup: see deploy/README.md).
# Needs the deploy key ~/.ssh/mzdesk_deploy on root@143.198.233.45.
set -euo pipefail
cd "$(dirname "$0")/.."
HOST=root@143.198.233.45
SSH="ssh -i $HOME/.ssh/mzdesk_deploy -o IdentitiesOnly=yes"
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
echo "== build frontend (base /mzdesk/)"
(cd web && MZDESK_BASE=/mzdesk/ npx vite build --outDir "$STAGE/static" --emptyOutDir --logLevel warn)
echo "== sync code -> /opt/mzdesk"
rsync -a --info=progress2 --human-readable --delete --exclude '__pycache__' --exclude 'static/' --exclude '.venv' \
  -e "$SSH" mzdesk requirements.txt pyproject.toml "$HOST:/opt/mzdesk/"
rsync -a --info=progress2 --human-readable --delete -e "$SSH" "$STAGE/static/" "$HOST:/opt/mzdesk/mzdesk/static/"
echo "== restart"
$SSH "$HOST" 'chown -R root:mzdesk /opt/mzdesk && chmod -R u=rwX,g=rX,o= /opt/mzdesk && systemctl restart mzdesk && sleep 4 && systemctl is-active mzdesk'
curl -s -o /dev/null -w "== live: https://bowenyang.org/mzdesk/api/spec -> %{http_code}\n" https://bowenyang.org/mzdesk/api/spec
