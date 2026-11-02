#!/bin/sh
# mzDesk installer for macOS and Linux (x86_64 and arm64).
#
#   curl -fsSL https://mzdesk.com/install.sh | sh
#
# Run it again to update. It never uses sudo and writes only to ~/.mzdesk, plus one link in ~/.local/bin or one
# line in your shell profile so that the `mzdesk` command is found.
#
# What it does, in order:
#   1. find the system          2. find uv (installs it into ~/.mzdesk/bin if there is none)
#   3. read version.json        4. download the wheel and check its sha256
#   5. make ~/.mzdesk/venv      6. install the wheel and its packages
#   7. write the launcher       8. put the launcher on your PATH
#
# Environment:
#   MZDESK_BASE_URL   where version.json and the wheel come from (default https://mzdesk.com)
#   MZDESK_VERSION    install this version instead of the newest one
#   MZDESK_UPDATING   set to 1 by `mzdesk update`: the last line then says "Updated to <version>"
set -eu

MZ="${HOME:-}/.mzdesk"
BASE="${MZDESK_BASE_URL:-https://mzdesk.com}"
BASE="${BASE%/}"
STEP="starting"          # named in the error message when something fails
TMP=""                   # scratch folder, removed on exit

say() { printf '%s\n' "$*"; }

# One clear sentence and the step that failed. The marker file tells the exit trap that the message is already out
# (fail may run inside $( ), where a plain variable would not survive).
fail() {
  printf 'mzDesk was not installed. Step "%s": %s\n' "$STEP" "$1" >&2
  [ -n "$TMP" ] && : > "$TMP/reported"
  exit 1
}

cleanup() {
  code=$?
  if [ "$code" -ne 0 ] && [ -n "$TMP" ] && [ ! -e "$TMP/reported" ]; then
    printf 'mzDesk was not installed. Step "%s" failed (the message above says why).\n' "$STEP" >&2
  fi
  [ -n "$TMP" ] && rm -rf "$TMP"
  return 0
}

# Value of a key in version.json. The file has one key per line, but this also reads it when it is on one line.
json_get() {
  printf '%s' "$JSON" | tr ',{}' '\n\n\n' |
    sed -n "s/^[[:space:]]*\"$1\"[[:space:]]*:[[:space:]]*\"\{0,1\}\([^\"]*\)\"\{0,1\}[[:space:]]*\$/\1/p" | head -n 1
}

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d ' ' -f 1
  elif command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -d ' ' -f 1
  elif command -v openssl >/dev/null 2>&1; then openssl dgst -sha256 "$1" | sed 's/^.*= *//'
  else fail "no tool to compute a sha256 checksum (sha256sum, shasum or openssl)."
  fi
}

# download URL FILE: a progress bar when a person is watching, quiet otherwise
download() {
  if [ -t 2 ]; then bar="-#"; else bar="-sS"; fi
  curl -fL "$bar" --retry 3 --connect-timeout 20 -o "$2" "$1"
}

main() {
  [ -n "${HOME:-}" ] || fail "the HOME variable is not set."
  command -v curl >/dev/null 2>&1 || fail "curl is not installed."
  mkdir -p "$MZ/cache" "$MZ/bin"
  TMP=$(mktemp -d "$MZ/cache/install.XXXXXX")

  # ---- 1. system
  STEP="finding your system"
  os=$(uname -s)
  arch=$(uname -m)
  case $os in
    Darwin) os=macos ;;
    Linux) os=linux ;;
    *) fail "this installer is for macOS and Linux; on Windows run: irm $BASE/install.ps1 | iex" ;;
  esac
  case $arch in
    x86_64 | amd64) arch=x86_64 ;;
    arm64 | aarch64) arch=arm64 ;;
    *) fail "the CPU type $arch is not supported (x86_64 and arm64 are)." ;;
  esac
  if [ "$os" = linux ] && (ldd --version 2>&1 || true) | grep -qi musl; then
    fail "Linux with musl (for example Alpine) is not supported; use a glibc system such as Ubuntu or Debian."
  fi
  say "System: $os $arch"

  # ---- 2. uv (it fetches Python and installs the packages)
  STEP="finding uv"
  if [ -x "$MZ/bin/uv" ]; then
    uv="$MZ/bin/uv"
  elif command -v uv >/dev/null 2>&1; then
    uv=$(command -v uv)
  else
    STEP="installing uv"
    say "Installing uv, the tool that sets up Python (astral.sh) ..."
    download https://astral.sh/uv/install.sh "$TMP/uv-install.sh" || fail "cannot download the uv installer from astral.sh (is the internet reachable?)."
    # UV_UNMANAGED_INSTALL: put uv into this folder and leave the PATH and the shell profiles alone
    UV_UNMANAGED_INSTALL="$MZ/bin" sh "$TMP/uv-install.sh" >"$TMP/uv.log" 2>&1 </dev/null || { cat "$TMP/uv.log" >&2; fail "the uv installer failed."; }
    uv="$MZ/bin/uv"
    [ -x "$uv" ] || fail "uv was not installed into $MZ/bin."
  fi

  # ---- 3. version.json
  STEP="reading version.json"
  JSON=$(curl -fsSL --retry 3 --connect-timeout 20 "$BASE/version.json") || fail "cannot download $BASE/version.json (is the internet reachable, and is the address right?)."
  latest=$(json_get version)
  if [ -n "${MZDESK_VERSION:-}" ] && [ "${MZDESK_VERSION}" != "$latest" ]; then
    # an older (or newer) version than the one version.json describes: its checksum sits next to the wheel
    version=$MZDESK_VERSION
    printf '%s' "$version" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$' || fail "MZDESK_VERSION must look like 0.5.4."
    wheel="mzdesk-$version-py3-none-any.whl"
    shafile=$(curl -fsSL --retry 3 --connect-timeout 20 "$BASE/download/$wheel.sha256") || fail "version $version is not available on $BASE."
    sha=$(printf '%s' "$shafile" | cut -d ' ' -f 1)
  else
    version=$latest
    wheel=$(json_get wheel)
    sha=$(json_get sha256)
  fi
  printf '%s' "$version" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$' || fail "version.json has no valid version."
  [ "$wheel" = "mzdesk-$version-py3-none-any.whl" ] || fail "version.json names an unexpected wheel file."
  printf '%s' "$sha" | grep -Eq '^[0-9a-f]{64}$' || fail "version.json has no valid sha256."

  # ---- 4. wheel
  STEP="downloading mzDesk $version"
  say "Downloading mzDesk $version ..."
  download "$BASE/download/$wheel" "$TMP/$wheel" || fail "cannot download the wheel from $BASE/download/$wheel."
  STEP="checking the download"
  got=$(sha256_of "$TMP/$wheel")
  [ "$got" = "$sha" ] || fail "the downloaded file does not match its sha256 checksum; try again, and report it if it repeats."

  # ---- 5. Python environment. A working one is kept, so an update only replaces the packages that changed.
  STEP="creating the Python environment"
  venv="$MZ/venv"
  if [ -x "$venv/bin/python" ] && "$venv/bin/python" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' >/dev/null 2>&1; then
    say "Using the Python environment in $venv"
  else
    say "Setting up Python 3.12 in $venv ..."
    rm -rf "$venv"
    # --managed-python: Python comes from uv, so the environment does not depend on a Python of the system
    "$uv" venv "$venv" --python 3.12 --managed-python </dev/null || fail "uv could not create the Python environment (it needs the internet to fetch Python)."
  fi

  # ---- 6. packages
  STEP="installing mzDesk and its packages"
  say "Installing mzDesk $version and its packages (about 1 GB; a few minutes) ..."
  if [ "$os" = linux ]; then
    # PyPI's PyTorch for Linux brings gigabytes of CUDA libraries; the CPU build is all mzDesk needs
    set -- --extra-index-url https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match
  else
    set --
  fi
  # --reinstall-package: a rebuilt wheel of the same version replaces the old one; the other packages stay if they still fit
  "$uv" pip install --python "$venv/bin/python" --reinstall-package mzdesk "$@" "$TMP/$wheel" </dev/null ||
    fail "uv could not install the packages (see the message above; the internet is needed)."

  # ---- 7. launcher
  STEP="writing the launcher"
  cat > "$MZ/bin/mzdesk" <<'EOF'
#!/bin/sh
# Written by the mzDesk installer: runs mzDesk from ~/.mzdesk/venv and opens it in the browser.
exec "$HOME/.mzdesk/venv/bin/python" -m mzdesk --open "$@"
EOF
  chmod 755 "$MZ/bin/mzdesk"

  # ---- 8. PATH
  STEP="putting mzdesk on your PATH"
  new_terminal=0
  local_bin="$HOME/.local/bin"
  case ":$PATH:" in
    *":$MZ/bin:"*) ;;                                  # already there (an update, or a shell that has read the profile line)
    *":$local_bin:"*)
      # ~/.local/bin is already on the PATH: a link is enough (nothing of the user's own is replaced)
      if [ -L "$local_bin/mzdesk" ] || [ ! -e "$local_bin/mzdesk" ]; then
        mkdir -p "$local_bin"
        ln -sf "$MZ/bin/mzdesk" "$local_bin/mzdesk"
      else
        add_to_profile
      fi ;;
    *)
      # otherwise, once: one line in the profile of the user's shell
      add_to_profile ;;
  esac
  [ "$new_terminal" = 0 ] || say "Open a new terminal window so that the mzdesk command is found."

  if [ "${MZDESK_UPDATING:-}" = 1 ]; then
    say "Updated to $version"
  else
    say "mzDesk $version is installed. Run: mzdesk"
  fi
}

# Add ~/.mzdesk/bin to the PATH in the profile of the user's login shell, unless a line for it is there already.
add_to_profile() {
  case ${SHELL:-} in
    */zsh) profile="${ZDOTDIR:-$HOME}/.zshrc" ;;
    */bash)
      if [ "$(uname -s)" != Darwin ]; then profile="$HOME/.bashrc"          # a login shell on macOS reads .bash_profile, then .profile
      elif [ -f "$HOME/.bash_profile" ] || [ ! -f "$HOME/.profile" ]; then profile="$HOME/.bash_profile"
      else profile="$HOME/.profile"
      fi ;;
    */fish) profile="$HOME/.config/fish/config.fish" ;;
    *) profile="$HOME/.profile" ;;
  esac
  if [ -f "$profile" ] && grep -Fq '.mzdesk/bin' "$profile"; then return 0; fi
  mkdir -p "$(dirname "$profile")"
  case $profile in
    */config.fish) line='set -gx PATH "$HOME/.mzdesk/bin" $PATH' ;;
    *) line='export PATH="$HOME/.mzdesk/bin:$PATH"' ;;
  esac
  printf '\n# mzDesk\n%s\n' "$line" >> "$profile"
  say "Added ~/.mzdesk/bin to the PATH in $profile"
  new_terminal=1
}

trap cleanup EXIT
# main runs last, so a download that stops half way never starts installing
main
