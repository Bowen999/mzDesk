#!/bin/bash
# Update the live site (https://mzdesk.com/) from this checkout: build the frontend, sync the code, install the service unit and
# the nginx files when they changed, restart, then publish the release for the local version (the wheel and version.json that
# `curl -fsSL https://mzdesk.com/install.sh | sh` installs; see "Releases and the install command" in deploy/README.md).
# Models and data under /var/lib/mzdesk are not touched.
# First-time setup and the move to mzdesk.com: deploy/README.md, deploy/cutover.sh. Needs the deploy key ~/.ssh/mzdesk_deploy.
set -euo pipefail
cd "$(dirname "$0")/.."
[ -z "${MZDESK_RELEASE_ALLOW_STUB_STATIC:-}" ] || { echo "!! MZDESK_RELEASE_ALLOW_STUB_STATIC is for tests of the installer; a deploy builds the real frontend."; exit 1; }
HOST=root@143.198.233.45
SSH="ssh -i $HOME/.ssh/mzdesk_deploy -o IdentitiesOnly=yes"
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
# The frontend is built for the root of mzdesk.com. Until the server has that certificate (deploy/cutover.sh), the live
# https://bowenyang.org/mzdesk/ still serves the old build, so deploying now would break it.
if ! $SSH "$HOST" 'test -f /etc/letsencrypt/live/mzdesk.com/fullchain.pem'; then
  echo "!! The server has no certificate for mzdesk.com yet: point the DNS records at ${HOST#*@}, then run deploy/cutover.sh (it deploys too)."
  exit 1
fi
echo "== build frontend"
(cd web && npx vite build --outDir "$STAGE/static" --emptyOutDir --logLevel warn)
echo "== build release (wheel + version.json -> dist/release)"
MZDESK_STATIC_DIR="$STAGE/static" deploy/build-release.sh
echo "== sync code -> /opt/mzdesk"
rsync -a --info=progress2 --human-readable --delete --exclude '__pycache__' --exclude 'static/' --exclude '.venv' \
  -e "$SSH" mzdesk requirements.txt pyproject.toml "$HOST:/opt/mzdesk/"
rsync -a --info=progress2 --human-readable --delete -e "$SSH" "$STAGE/static/" "$HOST:/opt/mzdesk/mzdesk/static/"
echo "== service unit, nginx, restart"
$SSH "$HOST" 'mkdir -p /root/mzdesk-deploy'
# the release waits in the staging folder; remote-install.sh publishes it to /var/lib/mzdesk/releases once the site answers
rsync -a --delete --info=progress2 --human-readable -e "$SSH" dist/release/ "$HOST:/root/mzdesk-deploy/release/"
rsync -a -e "$SSH" deploy/mzdesk.service deploy/00-mzdesk-limits.conf deploy/nginx-mzdesk.conf deploy/requirements-server.txt "$HOST:/root/mzdesk-deploy/"
$SSH "$HOST" 'bash -s' < deploy/remote-install.sh
if [ "$(dig +short mzdesk.com A | head -1)" = "${HOST#*@}" ]; then
  curl -s -o /dev/null -w "== live: https://mzdesk.com/api/health -> %{http_code}\n" https://mzdesk.com/api/health
  curl -s -o /dev/null -w "== live: https://mzdesk.com/install.sh -> %{http_code}\n" https://mzdesk.com/install.sh
  echo "== live: https://mzdesk.com/version.json"; curl -s https://mzdesk.com/version.json
else
  echo "== mzdesk.com does not point at ${HOST#*@} yet (see deploy/README.md); the app answers on the server only"
fi
