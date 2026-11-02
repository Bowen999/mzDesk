"""The installers (mzdesk/install/install.sh, install.ps1) and deploy/build-release.sh.
install.sh is run for real against a local server with a fake `uv`, in a throw-away HOME: no internet, no real install.
install.ps1 cannot be run here (no PowerShell): it gets syntax-free text checks for the known Windows PowerShell 5.1 pitfalls. Not verified on Windows."""
import hashlib
import http.server
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SH = ROOT / "mzdesk" / "install" / "install.sh"
PS1 = ROOT / "mzdesk" / "install" / "install.ps1"
DEPLOY = ROOT / "deploy"
needs_posix = pytest.mark.skipif(sys.platform == "win32" or not shutil.which("curl"), reason="install.sh needs a POSIX shell and curl")


# ---------------------------------------------------------------- the text of the scripts
@pytest.mark.parametrize("shell", ["sh", "bash", "dash"])
def test_install_sh_parses(shell):
    if not shutil.which(shell):
        pytest.skip(f"{shell} is not installed")
    subprocess.run([shell, "-n", str(SH)], check=True)


@pytest.mark.parametrize("script", ["build-release.sh", "deploy.sh", "remote-install.sh"])
def test_deploy_scripts_parse(script):
    if not shutil.which("bash"):
        pytest.skip("bash is not installed")
    subprocess.run(["bash", "-n", str(DEPLOY / script)], check=True)


def test_shellcheck_is_clean():
    if not shutil.which("shellcheck"):
        pytest.skip("shellcheck is not installed")
    subprocess.run(["shellcheck", "-s", "sh", str(SH)], check=True)
    subprocess.run(["shellcheck", "-s", "bash", str(DEPLOY / "build-release.sh"), str(DEPLOY / "deploy.sh")], check=True)


@pytest.mark.parametrize("path", [SH, PS1])
def test_scripts_are_plain_ascii_and_never_ask_for_root(path):
    text = path.read_bytes().decode("ascii")                    # raises on any other character: Windows PowerShell 5.1 reads a script without BOM as ANSI
    code = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))
    assert not re.search(r"\bsudo\b|\bsu\b|RunAs|-Verb", code)


@pytest.mark.parametrize("path", [SH, PS1])
def test_every_step_of_the_contract_is_in_the_script(path):
    text = path.read_text()
    for needle in ("MZDESK_BASE_URL", "MZDESK_VERSION", "MZDESK_UPDATING", "https://mzdesk.com", "astral.sh/uv/install", "UV_UNMANAGED_INSTALL",
                   "version.json", "/download/", "sha256", "--python", "3.12", "--managed-python", "pip", "install", "--reinstall-package", "mzdesk",
                   "-m mzdesk --open", "PATH", "is installed. Run: mzdesk", "Updated to", ".mzdesk"):
        assert needle in text, needle
    for step in ("finding your system", "finding uv", "installing uv", "reading version.json", "checking the download",
                 "creating the Python environment", "installing mzDesk and its packages", "writing the launcher", "on your PATH"):
        assert step in text, step


def test_install_sh_uses_the_cpu_pytorch_on_linux_only():
    text = SH.read_text()
    assert "https://download.pytorch.org/whl/cpu" in text and "--index-strategy unsafe-best-match" in text
    assert "https://download.pytorch.org/whl/cpu" not in PS1.read_text()          # PyPI's Windows PyTorch is the CPU build already


def test_install_ps1_avoids_the_windows_powershell_5_1_pitfalls():
    code = "\n".join(line for line in PS1.read_text().splitlines() if not line.lstrip().startswith("#"))
    assert "&&" not in code and "||" not in code and "??" not in code           # PowerShell 7 only
    assert "-AdditionalChildPath" not in code                                    # Join-Path with three parts: PowerShell 6+
    assert "$ProgressPreference = 'SilentlyContinue'" in code                    # Invoke-WebRequest is very slow with the progress bar
    assert "Tls12" in code                                                       # 5.1 may not offer TLS 1.2
    assert "$ErrorActionPreference = 'Stop'" in code
    assert "-UseBasicParsing" in code
    assert not re.search(r"^\s*exit\b", code, re.M)                              # through iex, `exit` would close the user's window
    assert "Set-ExecutionPolicy" not in code
    assert re.search(r"\[Environment\]::SetEnvironmentVariable\('Path', .*'User'\)", code)       # the user's PATH, no administrator rights
    assert not re.search(r'"[^"\n]*\$(?!env:|global:|script:)\w+:[^"\n]*"', code.replace("`\"", "~"))     # "$name: text" is read as a drive-qualified variable
    assert code.count("{") == code.count("}") and code.count("(") == code.count(")")


# ---------------------------------------------------------------- install.sh against a local server and a fake uv
class Releases:
    """Serves version.json at / and the wheels under /download/ from a folder, like mzdesk/api/release.py."""

    def __init__(self, folder: Path):
        folder_ = folder

        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                name = self.path.split("?")[0]
                name = name[len("/download/"):] if name.startswith("/download/") else name.lstrip("/")
                path = folder_ / name
                if "/" in name or not path.is_file():
                    self.send_error(404)
                    return
                body = path.read_bytes()
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def log_message(self, *args):
                pass

        self.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.server.server_address[1]}"
        threading.Thread(target=self.server.serve_forever, daemon=True).start()

    def close(self):
        self.server.shutdown()
        self.server.server_close()


def make_wheel(folder: Path, version: str) -> str:
    name = f"mzdesk-{version}-py3-none-any.whl"
    with zipfile.ZipFile(folder / name, "w") as z:
        z.writestr("mzdesk/__init__.py", f'__version__ = "{version}"\n')
    digest = hashlib.sha256((folder / name).read_bytes()).hexdigest()
    (folder / f"{name}.sha256").write_text(f"{digest}  {name}\n")
    return digest


FAKE_UV = """#!/bin/sh
echo "$*" >> "$FAKE_UV_LOG"
if [ "$1" = venv ]; then
  mkdir -p "$2/bin"
  printf '#!/bin/sh\\nexit 0\\n' > "$2/bin/python"
  chmod +x "$2/bin/python"
fi
"""


@pytest.fixture
def site(tmp_path):
    """A release folder (9.9.9 is the newest, 9.9.8 an older one), a server for it, a fake uv, and a throw-away HOME."""
    folder = tmp_path / "releases"
    folder.mkdir()
    sha = make_wheel(folder, "9.9.9")
    make_wheel(folder, "9.9.8")
    info = {"version": "9.9.9", "released": "2026-11-02", "wheel": "mzdesk-9.9.9-py3-none-any.whl", "sha256": sha, "size": 1, "python": ">=3.10",
            "notes": "https://mzdesk.com/docs/changelog"}
    (folder / "version.json").write_text(json.dumps(info, indent=2) + "\n")
    fakebin = tmp_path / "fakebin"
    fakebin.mkdir()
    (fakebin / "uv").write_text(FAKE_UV)
    (fakebin / "uv").chmod(0o755)
    home = tmp_path / "home"
    home.mkdir()
    server = Releases(folder)
    log = tmp_path / "uv.log"

    def run(*, base=None, path=None, **env):
        """Run install.sh the way `curl ... | sh` does, with nothing from the real environment but a minimal PATH."""
        full = {"HOME": str(home), "PATH": path or f"{fakebin}:/usr/bin:/bin", "SHELL": "/bin/zsh", "FAKE_UV_LOG": str(log),
                "MZDESK_BASE_URL": base or server.url, **env}
        return subprocess.run(["sh"], input=SH.read_bytes(), env=full, capture_output=True, timeout=120)

    yield type("Site", (), {"run": staticmethod(run), "home": home, "log": log, "folder": folder, "fakebin": fakebin, "sha": sha})
    server.close()


@needs_posix
def test_install_runs_every_step_and_prints_the_run_line(site):
    r = site.run()
    out, err = r.stdout.decode(), r.stderr.decode()
    assert r.returncode == 0, err
    assert out.rstrip().splitlines()[-1] == "mzDesk 9.9.9 is installed. Run: mzdesk"
    calls = site.log.read_text().splitlines()
    venv = site.home / ".mzdesk" / "venv"
    assert calls[0] == f"venv {venv} --python 3.12 --managed-python"
    install = calls[1].split()
    assert install[:5] == ["pip", "install", "--python", f"{venv}/bin/python", "--reinstall-package"] and install[5] == "mzdesk"
    assert install[-1].endswith("/mzdesk-9.9.9-py3-none-any.whl")
    if sys.platform.startswith("linux"):
        assert "--extra-index-url" in install and "--index-strategy" in install and "unsafe-best-match" in install
    else:
        assert "--extra-index-url" not in install
    launcher = site.home / ".mzdesk" / "bin" / "mzdesk"
    assert os.access(launcher, os.X_OK)
    assert launcher.read_text().splitlines()[-1] == 'exec "$HOME/.mzdesk/venv/bin/python" -m mzdesk --open "$@"'
    assert (site.home / ".zshrc").read_text().count(".mzdesk/bin") == 1                  # not on the PATH yet: one profile line
    assert "new terminal" in out
    assert list((site.home / ".mzdesk" / "cache").iterdir()) == []                        # the scratch folder is gone
    assert not (site.home / ".local").exists() and not (site.home / "tools").exists()     # nothing outside ~/.mzdesk and the profile


@needs_posix
def test_running_it_again_updates_and_changes_nothing_else(site):
    assert site.run().returncode == 0
    r = site.run(path=f"{site.home}/.local/bin:{site.fakebin}:/usr/bin:/bin", MZDESK_UPDATING="1")
    assert r.returncode == 0, r.stderr.decode()
    assert r.stdout.decode().rstrip().splitlines()[-1] == "Updated to 9.9.9"
    calls = site.log.read_text().splitlines()
    assert [c.split()[0] for c in calls] == ["venv", "pip", "pip"]                      # the working environment is kept, the wheel is installed again
    link = site.home / ".local" / "bin" / "mzdesk"
    assert link.is_symlink() and link.resolve() == (site.home / ".mzdesk" / "bin" / "mzdesk").resolve()
    assert (site.home / ".zshrc").read_text().count(".mzdesk/bin") == 1                  # the profile line is not added twice
    assert site.run(path=f"{site.home}/.local/bin:{site.fakebin}:/usr/bin:/bin").returncode == 0      # and a third time


@needs_posix
def test_a_pinned_version_uses_the_checksum_next_to_its_wheel(site):
    r = site.run(MZDESK_VERSION="9.9.8")
    assert r.returncode == 0, r.stderr.decode()
    assert r.stdout.decode().rstrip().splitlines()[-1] == "mzDesk 9.9.8 is installed. Run: mzdesk"
    assert site.log.read_text().splitlines()[1].endswith("/mzdesk-9.9.8-py3-none-any.whl")


@needs_posix
def test_a_wrong_checksum_installs_nothing(site):
    info = json.loads((site.folder / "version.json").read_text())
    (site.folder / "version.json").write_text(json.dumps({**info, "sha256": "0" * 64}))
    r = site.run()
    err = r.stderr.decode()
    assert r.returncode != 0
    assert 'Step "checking the download"' in err and "checksum" in err and err.count("\n") == 1      # one sentence, and the step
    assert not site.log.exists() and not (site.home / ".mzdesk" / "bin" / "mzdesk").exists()


@needs_posix
@pytest.mark.parametrize("env, step", [
    ({"MZDESK_BASE_URL": "http://127.0.0.1:1"}, "reading version.json"),                 # nothing listens there
    ({"MZDESK_VERSION": "9.9"}, "reading version.json"),
    ({"MZDESK_VERSION": "7.7.7"}, "reading version.json"),                                # not on the server
])
def test_failures_say_what_failed(site, env, step):
    r = site.run(**env)
    last = r.stderr.decode().strip().splitlines()[-1]                                   # curl may print its own line before ours
    assert r.returncode != 0 and last.startswith("mzDesk was not installed.") and f'Step "{step}"' in last
    assert not (site.home / ".mzdesk" / "bin" / "mzdesk").exists()


@needs_posix
def test_a_failing_uv_is_reported_with_its_step(site):
    (site.fakebin / "uv").write_text('#!/bin/sh\necho "no network" >&2\nexit 2\n')
    r = site.run()
    assert r.returncode != 0 and 'Step "creating the Python environment"' in r.stderr.decode()


# ---------------------------------------------------------------- deploy/build-release.sh
def test_build_release_makes_a_checked_wheel_and_version_json(tmp_path):
    if not (shutil.which("bash") and shutil.which("rsync")):
        pytest.skip("bash and rsync are needed")
    if not list((ROOT / "mzdesk" / "lib").glob("MetFragCommandLine-*.jar")):
        pytest.skip("the MetFrag jar is not in this checkout (it is not in git)")
    if not subprocess.run([sys.executable, "-c", "import build, setuptools, wheel"], capture_output=True).returncode == 0 and not shutil.which("uv"):
        pytest.skip("no wheel builder (pip install build, or uv)")
    out = tmp_path / "release"
    env = {**os.environ, "MZDESK_RELEASE_OUT": str(out), "MZDESK_RELEASE_ALLOW_STUB_STATIC": "1", "PYTHON": sys.executable}
    r = subprocess.run(["bash", str(DEPLOY / "build-release.sh")], env=env, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    info = json.loads((out / "version.json").read_text())
    wheel = out / info["wheel"]
    assert re.fullmatch(r"mzdesk-\d+\.\d+\.\d+-py3-none-any\.whl", info["wheel"]) and info["wheel"] == f"mzdesk-{info['version']}-py3-none-any.whl"
    assert info["sha256"] == hashlib.sha256(wheel.read_bytes()).hexdigest() and info["size"] == wheel.stat().st_size
    assert info["python"] == ">=3.10" and re.fullmatch(r"\d{4}-\d{2}-\d{2}", info["released"]) and info["notes"].endswith("/docs/changelog")
    assert (out / f"{info['wheel']}.sha256").read_text().split() == [info["sha256"], info["wheel"]]
    names = zipfile.ZipFile(wheel).namelist()
    for needed in ("mzdesk/static/index.html", "mzdesk/install/install.sh", "mzdesk/install/install.ps1", "mzdesk/api/release.py"):
        assert needed in names
    assert any(n.startswith("mzdesk/lib/MetFragCommandLine-") and n.endswith(".jar") for n in names)
    assert not any(n.split("/")[0] in ("benchmark", "tests", "web", "deploy") or "__pycache__" in n or n.endswith(".sqlite") for n in names)


def test_build_release_stops_without_the_frontend(tmp_path):
    if not (shutil.which("bash") and shutil.which("rsync")):
        pytest.skip("bash and rsync are needed")
    if not list((ROOT / "mzdesk" / "lib").glob("MetFragCommandLine-*.jar")):
        pytest.skip("the MetFrag jar is not in this checkout (it is not in git)")
    env = {**os.environ, "MZDESK_RELEASE_OUT": str(tmp_path / "release"), "MZDESK_STATIC_DIR": str(tmp_path / "no-such-build")}
    env.pop("MZDESK_RELEASE_ALLOW_STUB_STATIC", None)
    r = subprocess.run(["bash", str(DEPLOY / "build-release.sh")], env=env, capture_output=True, text=True, timeout=60)
    assert r.returncode != 0 and "frontend is not built" in r.stderr and "npm run build" in r.stderr
    assert not (tmp_path / "release").exists()
