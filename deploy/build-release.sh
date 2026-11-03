#!/bin/bash
# Build the release that https://mzdesk.com hands out to the local version, into dist/release/:
#   mzdesk-X.Y.Z-py3-none-any.whl   the wheel (X.Y.Z from mzdesk/__init__.py)
#   mzdesk-X.Y.Z-py3-none-any.whl.sha256
#   version.json                    what install.sh and `mzdesk update` read
# Nothing is uploaded and no server is contacted; deploy/deploy.sh uploads dist/release/. Safe to run again.
#
# The wheel is built from a copy of the tree in a temporary folder, so the checkout gets no build/ or egg-info folders.
# It needs the built frontend (cd web && npm run build, which fills mzdesk/static) and the MetFrag jar (mzdesk/lib/, not in git).
#
# Environment (all optional):
#   MZDESK_STATIC_DIR=<dir>              use this built frontend instead of mzdesk/static (deploy.sh builds one into a temp folder)
#   MZDESK_RELEASE_OUT=<dir>             write there instead of dist/release
#   MZDESK_RELEASE_ALLOW_STUB_STATIC=1   FOR TESTS OF THE INSTALLER ONLY: when there is no built frontend (or no docs figures), build anyway
#                                        with a stub index.html in the copy and print a loud warning. Never deploy such a wheel.
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON=${PYTHON:-python3}
OUT=${MZDESK_RELEASE_OUT:-dist/release}
STATIC_SRC=${MZDESK_STATIC_DIR:-mzdesk/static}
STUB=${MZDESK_RELEASE_ALLOW_STUB_STATIC:-0}
die() { echo "!! $*" >&2; exit 1; }

VERSION=$(sed -n 's/^__version__ = "\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\)".*/\1/p' mzdesk/__init__.py)
[ -n "$VERSION" ] || die "cannot read __version__ (X.Y.Z) from mzdesk/__init__.py"
WHEEL="mzdesk-$VERSION-py3-none-any.whl"
echo "== mzDesk $VERSION"

# ---- what the wheel is built from
JAR=$(ls mzdesk/lib/MetFragCommandLine-*.jar 2>/dev/null | head -1 || true)
[ -n "$JAR" ] || die "the MetFrag jar is missing: mzdesk/lib/MetFragCommandLine-*.jar (it is not in git; copy it from a previous checkout or from the server's /opt/mzdesk/mzdesk/lib/)"
STUBBED=0
if [ ! -f "$STATIC_SRC/index.html" ]; then
  [ "$STUB" = 1 ] || die "the frontend is not built: $STATIC_SRC/index.html is missing. Run: cd web && npm run build"
  echo "!! WARNING: no built frontend. MZDESK_RELEASE_ALLOW_STUB_STATIC=1 puts a STUB page into the wheel."
  echo "!! WARNING: this wheel is for testing the installer only. Do not deploy it."
  STUBBED=1
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/src" "$WORK/out"
# the package without caches and without the frontend, then the frontend from where it was built
rsync -a --exclude '__pycache__' --exclude '*.pyc' --exclude '.DS_Store' --exclude '/mzdesk/static/' \
  mzdesk pyproject.toml requirements.txt LICENSE "$WORK/src/"
if [ "$STUBBED" = 0 ]; then
  rsync -a --exclude '.DS_Store' "$STATIC_SRC"/ "$WORK/src/mzdesk/static/"
else
  mkdir -p "$WORK/src/mzdesk/static"
  echo '<!doctype html><title>mzDesk stub</title><p>STUB frontend: this wheel is only for testing the installer.</p>' > "$WORK/src/mzdesk/static/index.html"
fi

# ---- build. --no-isolation needs no network (setuptools is already installed); uv build fetches setuptools from PyPI if it has none cached.
echo "== build the wheel"
if "$PYTHON" -c 'import build, setuptools, wheel' 2>/dev/null; then
  "$PYTHON" -m build --wheel --no-isolation --outdir "$WORK/out" "$WORK/src" >"$WORK/build.log" 2>&1 || { tail -30 "$WORK/build.log" >&2; die "the wheel build failed"; }
elif command -v uv >/dev/null 2>&1; then
  uv build --wheel --out-dir "$WORK/out" "$WORK/src" >"$WORK/build.log" 2>&1 || { tail -30 "$WORK/build.log" >&2; die "the wheel build failed"; }
else
  die "no wheel builder: install uv (https://docs.astral.sh/uv/), or run this script in a virtual environment that has: pip install build"
fi
[ "$(ls "$WORK/out" | wc -l | tr -d ' ')" = 1 ] && [ -f "$WORK/out/$WHEEL" ] || die "expected the wheel $WHEEL, got: $(ls "$WORK/out" | tr '\n' ' ')"

# ---- look inside: what must be there, what must not, what is big
echo "== check the wheel"
if ! "$PYTHON" - "$WORK/out/$WHEEL" "$VERSION" "$STUB" <<'PY'
import fnmatch, sys, zipfile

path, version, stub = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
z = zipfile.ZipFile(path)
files = {i.filename: i for i in z.infolist() if not i.filename.endswith("/")}
names = sorted(files)
errors, warnings = [], []

def has(pattern):
    return any(fnmatch.fnmatchcase(n, pattern) for n in names)

# must be there
need = ["mzdesk/__init__.py", "mzdesk/api/app.py", "mzdesk/api/release.py", "mzdesk/static/index.html", "mzdesk/lib/MetFragCommandLine-*.jar",
        "mzdesk/install/install.sh", "mzdesk/install/install.ps1", "mzdesk/docs/*.md", "mzdesk/docs/img/*.svg", "mzdesk/data/*",
        "mzdesk/lib/msfinder/mzdesk_msfinder.ini", "mzdesk/lib/msfinder/ThreadStaticPatch.cs"]
for pattern in need:
    if not has(pattern):
        (warnings if stub and pattern in ("mzdesk/docs/img/*.svg",) else errors).append(f"missing in the wheel: {pattern}")
meta = next((n for n in names if n.endswith(".dist-info/METADATA")), None)
text = z.read(meta).decode("utf-8") if meta else ""
if f"\nVersion: {version}\n" not in text:
    errors.append(f"METADATA does not say Version: {version}")
if "\nRequires-Dist: jdk4py" not in text:
    errors.append("METADATA does not require jdk4py (requirements.txt)")

# must not be there: anything outside the package, development folders, secrets, caches, databases
bad_parts = {"benchmark", "tests", "node_modules", ".git", "__pycache__", "paper", "research_inputs"}
bad_names = [".env", ".env.*", "*.pem", "*.key", "id_rsa*", "*.sqlite", "*.sqlite-*", "*.db", "*.pyc", ".DS_Store", "*.log", "*secret*", "*password*", "*.egg-info"]
for n in names:
    parts = n.split("/")
    if parts[0] != "mzdesk" and not parts[0].endswith(".dist-info"):
        errors.append(f"outside the package: {n}")
    elif bad_parts & set(parts) or any(fnmatch.fnmatchcase(parts[-1].lower(), p) for p in bad_names):
        errors.append(f"does not belong in a release: {n}")

# worth a look: the big files, and data whose licence is not settled
total = sum(i.file_size for i in files.values())
big = sorted(((i.file_size, n) for n, i in files.items() if i.file_size > 3_000_000), reverse=True)
for size, n in big:
    if not fnmatch.fnmatchcase(n, "mzdesk/lib/MetFragCommandLine-*.jar"):
        warnings.append(f"large file ({size / 1e6:.1f} MB): {n}")
if has("mzdesk/gnps/redu_cache/*"):
    warnings.append("the Pan-ReDU snapshot mzdesk/gnps/redu_cache/ is in the wheel; .gitignore says its redistribution licence is not confirmed")

groups = {}
for n, i in files.items():
    key = "/".join(n.split("/")[:2]) if n.startswith("mzdesk/") else n.split("/")[0]
    groups.setdefault(key, [0, 0])
    groups[key][0] += 1
    groups[key][1] += i.file_size
print(f"   {len(files)} files, {total / 1e6:.1f} MB unpacked")
for key, (count, size) in sorted(groups.items(), key=lambda kv: -kv[1][1])[:8]:
    print(f"   {size / 1e6:8.1f} MB  {count:5d} files  {key}")
for w in warnings:
    print(f"   NOTE: {w}")
for e in errors:
    print(f"   !! {e}", file=sys.stderr)
sys.exit(1 if errors else 0)
PY
then
  die "the wheel is not right (see above); nothing was written to $OUT"
fi

# ---- write the release: wheel, checksum, then version.json (from a temp name, so it never exists half written)
sha256_of() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d ' ' -f 1; else shasum -a 256 "$1" | cut -d ' ' -f 1; fi; }
SHA=$(sha256_of "$WORK/out/$WHEEL")
SIZE=$(wc -c < "$WORK/out/$WHEEL" | tr -d ' ')
mkdir -p "$OUT"
rm -f "$OUT"/mzdesk-*-py3-none-any.whl "$OUT"/mzdesk-*-py3-none-any.whl.sha256 "$OUT/version.json"   # an earlier release is not part of this one
cp "$WORK/out/$WHEEL" "$OUT/$WHEEL"
echo "$SHA  $WHEEL" > "$OUT/$WHEEL.sha256"
cat > "$OUT/version.json.new" <<EOF
{
  "version": "$VERSION",
  "released": "$(date -u +%Y-%m-%d)",
  "wheel": "$WHEEL",
  "sha256": "$SHA",
  "size": $SIZE,
  "python": ">=3.10",
  "notes": "https://mzdesk.com/docs/changelog"
}
EOF
mv "$OUT/version.json.new" "$OUT/version.json"
echo "== release in $OUT"
echo "   $WHEEL  $((SIZE / 1000000)) MB  sha256 $SHA"
[ "$STUBBED" = 0 ] || echo "!! WARNING: this wheel has a STUB frontend (MZDESK_RELEASE_ALLOW_STUB_STATIC=1): for installer tests only. Do not deploy."
