"""The `mzdesk` command: arguments (also behind the launcher's --open), the port search, `setup`, `update`, `version`. Nothing is installed
and no server is started; uvicorn, the browser and the installer are replaced. HOME is a temporary folder."""
import dataclasses
import http.server
import os
import socket
import subprocess
import sys
import threading

import pytest

from mzdesk import __version__, cli, update
from mzdesk.local import components


@pytest.fixture(autouse=True)
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    for key in list(os.environ):
        if key.startswith("MZDESK_"):
            monkeypatch.delenv(key)
    monkeypatch.setenv("MZDESK_NO_UPDATE_CHECK", "1")
    from mzdesk.local import prereq
    monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True: "ready")
    return tmp_path


# ---------------------------------------------------------------- arguments
@pytest.fixture
def started(monkeypatch):
    calls = []
    monkeypatch.setattr(cli, "run_server", lambda args: calls.append(args) or 0)
    return calls


def test_run_arguments(started):
    assert cli.main([]) == 0
    a = started[-1]
    assert (a.port, a.host, a.open, a.reload) == (None, "127.0.0.1", False, False)
    cli.main(["--port", "8599", "--host", "0.0.0.0", "--open", "--reload"])
    a = started[-1]
    assert (a.port, a.host, a.open, a.reload) == (8599, "0.0.0.0", True, True)


def test_bad_arguments_exit_with_usage(started, capsys):
    for argv in (["--port", "x"], ["--nope"], ["bogus"]):
        with pytest.raises(SystemExit) as exc:
            cli.main(argv)
        assert exc.value.code == 2
    assert not started


def test_version(capsys):
    assert cli.main(["version"]) == 0 and cli.main(["--version"]) == 0 and cli.main(["-V"]) == 0
    assert capsys.readouterr().out.splitlines() == [f"mzDesk {__version__}"] * 3


@pytest.mark.parametrize("argv,needle", [(["--open", "version"], f"mzDesk {__version__}"), (["--open", "setup", "--list"], "Components (size is the download):"),
                                        (["--open", "--open", "setup", "--check", "ccs"], "ccs"), (["--open", "update", "--check"], "latest 99.0.0")])
def test_the_launchers_open_flag_does_not_hide_the_subcommand(argv, needle, started, monkeypatch, capsys):
    monkeypatch.setattr(update, "refresh", lambda timeout=3: {"latest": "99.0.0"})
    cli.main(argv)
    assert needle in capsys.readouterr().out and not started                              # the subcommand ran, the server did not start


def test_open_still_means_open_the_browser_for_the_server(started):
    cli.main(["--open"])
    assert started[-1].open is True
    cli.main(["--open", "--port", "9000"])
    assert started[-1].port == 9000


# ---------------------------------------------------------------- ports
def free_pair():
    """Two ports in a row that are free now (the search walks up from the first)."""
    for _ in range(50):
        with socket.socket() as s:
            s.bind(("127.0.0.1", 0))
            port = s.getsockname()[1]
        if port < 65000 and cli.port_is_free("127.0.0.1", port) and cli.port_is_free("127.0.0.1", port + 1):
            return port
    pytest.skip("no free pair of ports")


def test_port_search_skips_busy_ports():
    port = free_pair()
    assert cli.find_port("127.0.0.1", port, port + 1) == port
    with socket.socket() as busy:
        busy.bind(("127.0.0.1", port))
        busy.listen()
        assert not cli.port_is_free("127.0.0.1", port)
        assert cli.find_port("127.0.0.1", port, port + 1) == port + 1
        with socket.socket() as busy2:
            busy2.bind(("127.0.0.1", port + 1))
            busy2.listen()
            assert cli.find_port("127.0.0.1", port, port + 1) is None


def test_default_port_range():
    assert (cli.DEFAULT_PORT, cli.LAST_PORT) == (8731, 8750)


def run_args(**kw):
    ns = dict(port=None, host="127.0.0.1", open=False, reload=False)
    ns.update(kw)
    return type("Args", (), ns)()


@pytest.fixture
def fake_uvicorn(monkeypatch):
    calls = []
    import uvicorn
    monkeypatch.setattr(uvicorn, "run", lambda *a, **k: calls.append((a, k)))
    return calls


def test_server_starts_on_the_first_free_port_and_prints_the_address(fake_uvicorn, monkeypatch, capsys):
    monkeypatch.setattr(cli, "find_port", lambda host: 8734)
    assert cli.run_server(run_args()) == 0
    (args, kw) = fake_uvicorn[0]
    assert args == ("mzdesk.api.app:app",) and kw == {"host": "127.0.0.1", "port": 8734, "reload": False}
    assert f"mzDesk {__version__} is running at http://127.0.0.1:8734" in capsys.readouterr().out


def test_an_explicit_port_is_used_as_given(fake_uvicorn, monkeypatch):
    monkeypatch.setattr(cli, "find_port", lambda host: pytest.fail("no search when --port is given"))
    cli.run_server(run_args(port=8599, reload=True))
    assert fake_uvicorn[0][1] == {"host": "127.0.0.1", "port": 8599, "reload": True}


def test_all_ports_busy_is_a_clear_error(fake_uvicorn, monkeypatch, capsys):
    monkeypatch.setattr(cli, "find_port", lambda host: None)
    assert cli.run_server(run_args()) == 1 and not fake_uvicorn
    assert "8731-8750 are all in use" in capsys.readouterr().err


def test_the_start_prints_the_update_line_from_the_cache(fake_uvicorn, monkeypatch, capsys):
    monkeypatch.setattr(cli, "find_port", lambda host: 8731)
    monkeypatch.setattr(update, "terminal_line", lambda: "mzDesk 99.0.0 is available (you have 0.5.3). Update: mzdesk update")
    monkeypatch.setattr(update, "check_in_background", lambda: None)
    cli.run_server(run_args())
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "mzDesk 99.0.0 is available (you have 0.5.3). Update: mzdesk update" and "is running at" in out[1]


def test_open_opens_the_browser_once_the_server_answers(monkeypatch):
    opened = []
    monkeypatch.setattr(cli.webbrowser, "open", lambda url: opened.append(url))
    with socket.socket() as srv:
        srv.bind(("127.0.0.1", 0))
        srv.listen()
        port = srv.getsockname()[1]
        cli._open_when_up("127.0.0.1", port, wait_s=5)
    assert opened == [f"http://127.0.0.1:{port}"]
    opened.clear()
    cli._open_when_up("127.0.0.1", free_pair(), wait_s=0.5)                               # nothing listens: no browser
    assert opened == []
    assert cli.browser_url("0.0.0.0", 8731) == "http://localhost:8731"


def test_open_flag_starts_the_opener_thread(fake_uvicorn, monkeypatch):
    seen = []
    monkeypatch.setattr(cli, "find_port", lambda host: 8732)
    monkeypatch.setattr(cli, "_open_when_up", lambda host, port: seen.append((host, port)))
    cli.run_server(run_args(open=True))
    deadline = 50
    while not seen and deadline:
        threading.Event().wait(0.05)
        deadline -= 1
    assert seen == [("127.0.0.1", 8732)]
    seen.clear()
    cli.run_server(run_args())
    threading.Event().wait(0.2)
    assert seen == []


# ---------------------------------------------------------------- setup
def fake_registry(monkeypatch, installed=(), ran=None):
    ran = [] if ran is None else ran
    state = set(installed)
    reg = {}
    for cid, comp in components.REGISTRY.items():
        reg[cid] = dataclasses.replace(comp, check=lambda cid=cid: cid in state,
                                       run=lambda ctx, cid=cid: (ctx.log(f"installing {cid}"), ran.append(cid), state.add(cid)))
    monkeypatch.setattr(components, "REGISTRY", reg)
    return ran


def test_setup_list_and_check(monkeypatch, capsys):
    fake_registry(monkeypatch, installed=["molgpka"])
    assert cli.main(["setup", "--list"]) == 0
    out = capsys.readouterr().out
    assert "molgpka" in out and "installed" in out and "2.3 GB" in out and "Spectral libraries" in out
    assert [line.split()[0] for line in out.splitlines()[1:]] == list(components.REGISTRY)
    assert cli.main(["setup", "--check"]) == 1                                              # most are missing
    assert cli.main(["setup", "--check", "molgpka"]) == 0
    assert cli.main(["setup", "--check", "molgpka", "ccs"]) == 1
    capsys.readouterr()
    assert cli.main(["setup", "--list", "nope"]) == 2 and "Unknown component: nope" in capsys.readouterr().out


def test_setup_asks_per_component(monkeypatch, capsys):
    ran = fake_registry(monkeypatch)
    answers = iter(["y", "n", "", "YES"])
    asked = []
    args = cli.build_parsers()[1].parse_args(["ccs", "molgpka", "repo", "libraries"])
    code = cli.run_setup(args, ask=lambda q: asked.append(q) or next(answers))
    assert code == 0 and ran == ["libraries", "ccs"] or sorted(ran) == ["ccs", "libraries"]
    assert len(asked) == 4 and "(2.3 GB download)" in asked[3] and "Install CCS reference data" in asked[0]
    out = capsys.readouterr().out
    assert "molgpka: skipped." in out and "repo: skipped." in out and "installing ccs" in out


def test_setup_yes_installs_without_asking_and_streams_the_log(monkeypatch, capsys):
    ran = fake_registry(monkeypatch)
    assert cli.main(["setup", "--yes", "ccs", "molgpka"]) == 0
    assert ran == ["molgpka", "ccs"]                                                        # registry order
    out = capsys.readouterr().out
    assert "installing molgpka" in out and "molgpka: done" in out and "ccs: done" in out


def test_setup_all_with_yes_and_already_installed(monkeypatch, capsys):
    ran = fake_registry(monkeypatch, installed=["ccs"])
    assert cli.main(["setup", "--yes", "all"]) == 0
    assert "ccs" not in ran and set(ran) == set(components.REGISTRY) - {"ccs"} and "ccs: already installed." in capsys.readouterr().out


def test_setup_yes_without_names_does_nothing(monkeypatch, capsys):
    ran = fake_registry(monkeypatch)
    assert cli.main(["setup", "--yes"]) == 2 and ran == [] and "mzdesk setup --yes all" in capsys.readouterr().out


def test_setup_says_why_a_component_is_unavailable(monkeypatch, capsys):
    ran = fake_registry(monkeypatch)
    components.REGISTRY["msfinder"] = dataclasses.replace(components.REGISTRY["msfinder"], platforms=components.POSIX_ONLY)
    monkeypatch.setattr(components, "current_platform", lambda: "win32")
    assert cli.main(["setup", "--yes", "msfinder", "ccs"]) == 0
    out = capsys.readouterr().out
    assert "msfinder: not available on this computer (macOS / Linux only for now)." in out and ran == ["ccs"]
    cli.main(["setup", "--list"])
    assert "not available: macOS / Linux only for now" in capsys.readouterr().out


def test_setup_failure_exits_1_with_the_message(monkeypatch, capsys):
    reg = {cid: dataclasses.replace(c, check=lambda: False, run=lambda ctx: ctx.run([sys.executable, "-c", "print('no network'); raise SystemExit(4)"], what="the download"))
           for cid, c in components.REGISTRY.items()}
    monkeypatch.setattr(components, "REGISTRY", reg)
    assert cli.main(["setup", "--yes", "ccs"]) == 1
    out = capsys.readouterr().out
    assert "no network" in out and "the download failed (exit code 4)" in out and "ccs: error" in out


# ---------------------------------------------------------------- update
@pytest.fixture
def newest(monkeypatch):
    got = {"latest": "99.0.0"}
    monkeypatch.setattr(update, "refresh", lambda timeout=3: got)
    return got


def test_update_check_prints_current_and_latest(newest, capsys, monkeypatch):
    monkeypatch.setattr(cli, "run_installer", lambda base, out=print: pytest.fail("--check must not install"))
    assert cli.main(["update", "--check"]) == 0
    assert capsys.readouterr().out.strip() == f"Current {__version__}, latest 99.0.0. Update: mzdesk update"
    newest["latest"] = __version__
    cli.main(["update", "--check"])
    assert "up to date" in capsys.readouterr().out


def test_update_cannot_read_the_version(newest, capsys):
    newest["latest"] = None
    assert cli.main(["update"]) == 1 and "Could not read the newest version" in capsys.readouterr().out


def test_update_runs_the_installer_of_the_server(newest, monkeypatch):
    seen = []
    monkeypatch.setenv("MZDESK_BASE_URL", "http://localhost:9999/")
    monkeypatch.setattr(cli, "run_installer", lambda base, out=print: seen.append(base) or 0)
    assert cli.main(["update"]) == 0 and seen == ["http://localhost:9999"]
    newest["latest"] = __version__
    seen.clear()
    assert cli.main(["update"]) == 0 and seen == []


@pytest.fixture
def script_server():
    served = {}

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            served.setdefault("paths", []).append(self.path)
            body = served.get(self.path)
            self.send_response(200 if body else 404)
            self.end_headers()
            self.wfile.write(body or b"")

        def log_message(self, *a):
            pass

    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    served["base"] = f"http://127.0.0.1:{srv.server_port}"
    yield served
    srv.shutdown()


@pytest.mark.skipif(sys.platform == "win32", reason="runs install.sh")
def test_the_installer_runs_with_mzdesk_updating(script_server, home, capsys):
    out_file = home / "result.txt"
    script_server["/install.sh"] = f'echo "updating=$MZDESK_UPDATING" > "{out_file}"\nexit 7\n'.encode()
    assert cli.run_installer(script_server["base"]) == 7                                    # the script's exit code is ours
    assert out_file.read_text().strip() == "updating=1"
    assert f"Updating mzDesk with {script_server['base']}/install.sh" in capsys.readouterr().out


def test_the_installer_download_failure_is_reported(script_server, capsys):
    assert cli.run_installer(script_server["base"]) == 1                                     # 404
    assert "Could not download" in capsys.readouterr().out


def test_windows_starts_the_installer_detached_and_exits(script_server, monkeypatch, capsys):
    script_server["/install.ps1"] = b"Write-Host hi\n"
    popen = []
    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.setattr(subprocess, "Popen", lambda cmd, **kw: popen.append((cmd, kw)))
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: pytest.fail("Windows must not wait for the installer"))
    assert cli.run_installer(script_server["base"]) == 0
    cmd, kw = popen[0]
    assert cmd[:6] == ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", cmd[5]] and cmd[5].endswith("install.ps1")
    assert kw["env"]["MZDESK_UPDATING"] == "1" and kw["creationflags"] & 0x10 and kw["creationflags"] & 0x200       # new console, new process group
    assert capsys.readouterr().out.strip().splitlines()[-1] == "Updating in a new window; mzDesk will be closed."
    assert script_server["paths"] == ["/install.ps1"]
