#!/bin/bash
# Install, on the server, what the local version brings with it and the pages need: the spectral libraries (MassBank + GNPS), the CCS reference
# data, Mono + MS-FINDER, the FRIGID environment (CPU PyTorch) with its checkpoints, and the CFM-ID Docker image. Each step skips what is already
# in place, so it is safe to run again. The Python packages of the main venv come with deploy.sh (requirements-server.txt).
# About 7 GB downloaded, 10 GB on disk; FRIGID's checkpoint (2.6 GB, Zenodo) is the slow one. Needs the deploy key ~/.ssh/mzdesk_deploy.
set -euo pipefail
cd "$(dirname "$0")/.."
HOST=root@143.198.233.45
SSH="ssh -i $HOME/.ssh/mzdesk_deploy -o IdentitiesOnly=yes"
rsync -a -e "$SSH" deploy/requirements-server.txt "$HOST:/root/mzdesk-deploy/"
$SSH "$HOST" 'bash -s' <<'REMOTE'
set -euo pipefail
D=/var/lib/mzdesk
PY=/opt/mzdesk/.venv/bin/python
install -d -o mzdesk -g mzdesk -m 750 $D/data/libraries $D/data/ccs $D/tools $D/setup-logs
install -d -o mzdesk -g mzdesk -m 700 $D/cache/cfmid-work

# as_app <memory cap> <environment> <command...>: run as the service user, capped, with the output on screen
as_app() {
  local mem=$1 env=$2; shift 2
  systemd-run --pipe --wait --collect --quiet --uid=mzdesk --gid=mzdesk -p MemoryMax=$mem -p Nice=10 -p WorkingDirectory=/opt/mzdesk \
    -p Environment="HOME=$D PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 $env" "$@"
}

echo "== spectral libraries (MassBank 137 MB, GNPS 2.2 GB download; about 1.3 GB of index)"
if [ -d $D/data/libraries/MassBank ] && [ -d $D/data/libraries/GNPS ]; then echo "   there"; else
  as_app 2200M "MZDESK_LIBRARY_DIR=$D/data/libraries" $PY -m mzdesk.ms2.library_setup MassBank GNPS --delete-raw
fi

echo "== CCS reference data and model"
if [ -n "$(ls -A $D/data/ccs 2>/dev/null)" ]; then echo "   there"; else
  as_app 1200M "MZDESK_CCS_DIR=$D/data/ccs" $PY -m mzdesk.orthogonal.ccs_setup
fi

echo "== Mono (conda-forge: Ubuntu's package has no Roslyn csc) and MS-FINDER"
if [ ! -x $D/tools/mono/bin/csc ]; then
  /root/miniconda3/bin/conda create -y -q -p $D/tools/mono --override-channels -c conda-forge mono=6.12.0.199
  chown -R mzdesk:mzdesk $D/tools/mono
fi
as_app 1200M "MZDESK_MONO=$D/tools/mono/bin/mono MZDESK_MSFINDER_DIR=$D/tools/msfinder/MSFINDER_v3.73" $PY -m mzdesk.ms2.msfinder_setup

echo "== FRIGID: Python 3.10 environment with CPU PyTorch, code, checkpoints (2.9 GB from Zenodo, about an hour)"
if [ ! -x $D/tools/frigid/env/bin/python ]; then
  install -d -m 755 $D/tools/frigid
  /root/miniconda3/bin/conda create -y -q -p $D/tools/frigid/env --override-channels -c conda-forge python=3.10 antlr-python-runtime=4.9.3
  $D/tools/frigid/env/bin/python -m pip install -q --disable-pip-version-check torch==2.6.0 --index-url https://download.pytorch.org/whl/cpu
fi
chown -R mzdesk:mzdesk $D/tools/frigid
as_app 2G "MZDESK_FRIGID_DIR=$D/tools/frigid" $PY -m mzdesk.ms2.frigid_setup

echo "== CFM-ID Docker image"
IMAGE=wishartlab/cfmid:msml4.4.10_msrb1.1.15
docker image inspect $IMAGE >/dev/null 2>&1 && echo "   there" || docker pull --platform linux/amd64 $IMAGE
echo "== done"
REMOTE
