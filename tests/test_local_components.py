"""Local version: the component registry, installed() on fake folders, needs() per page, the install job (success, failure, cancel), the commit
archive download and the engine folder rules. No network (one local HTTP server), nothing is installed, HOME is a temporary folder."""
import dataclasses
import http.server
import io
import os
import sys
import tarfile
import threading
import time
from pathlib import Path

import pytest

from mzdesk.local import archive, components, jobs, paths, prereq, targets
from mzdesk.local.runner import Cancelled, Ctx, InstallError

IDS = ["libraries", "msfinder", "frigid", "cfmid", "neims", "molgpka", "ccs", "repo"]


@pytest.fixture(autouse=True)
def home(tmp_path, monkeypatch):
    """A fresh HOME for every test; no MZDESK_* folder switch leaks in from the developer's shell."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    for key in list(os.environ):
        if key.startswith("MZDESK_") and key != "MZDESK_HOSTED":
            monkeypatch.delenv(key)
    monkeypatch.delenv("MZDESK_HOSTED", raising=False)
    monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True: "ready")
    return tmp_path


def with_installed(monkeypatch, **flags):
    """The registry with installed() fixed per id (all others: not installed)."""
    reg = {cid: dataclasses.replace(c, check=lambda v=flags.get(cid, False): v) for cid, c in components.REGISTRY.items()}
    monkeypatch.setattr(components, "REGISTRY", reg)


# ---------------------------------------------------------------- registry
def test_registry_has_the_eight_components_in_install_order():
    assert list(components.REGISTRY) == IDS
    for comp in components.REGISTRY.values():
        assert comp.title and comp.size and callable(comp.check) and callable(comp.run) and set(comp.platforms) <= {"darwin", "linux", "win32"}
        assert comp.info().keys() == {"id", "title", "size", "installed", "supported", "ask", "note"}


def test_ask_flags_and_platforms_follow_the_contract():
    reg = components.REGISTRY
    assert [cid for cid in IDS if reg[cid].ask] == ["libraries", "msfinder", "frigid", "cfmid", "neims"]
    assert [cid for cid in IDS if "win32" in reg[cid].platforms] == ["libraries", "molgpka", "ccs"]
    assert all(set(reg[cid].platforms) == {"darwin", "linux"} for cid in ("msfinder", "frigid", "cfmid", "neims", "repo"))


def test_unsupported_platform_is_reported_with_the_reason(monkeypatch):
    monkeypatch.setattr(components, "current_platform", lambda: "win32")
    info = {i["id"]: i for i in components.status()}
    assert info["msfinder"]["supported"] is False and info["msfinder"]["installed"] is False and info["msfinder"]["note"] == "macOS / Linux only for now"
    assert info["libraries"]["supported"] is True
    with pytest.raises(InstallError, match="macOS / Linux only"):
        components.plan(["msfinder"])


def test_plan_orders_by_the_registry_and_adds_molgpka_before_repo(monkeypatch):
    with_installed(monkeypatch)
    assert components.plan(["repo", "libraries"]) == ["libraries", "molgpka", "repo"]
    with_installed(monkeypatch, molgpka=True)
    assert components.plan(["repo"]) == ["repo"]
    with pytest.raises(InstallError, match="Unknown component"):
        components.plan(["nope"])


# ---------------------------------------------------------------- installed() on fake folders
def touch(path, size=0, text=""):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    if size:
        os.truncate(path, size)                 # a sparse file: the engines only compare the size
    return path


def test_simple_components_follow_their_folders(home, monkeypatch):
    reg = components.REGISTRY
    assert not any(reg[c].installed() for c in ("libraries", "neims", "molgpka", "ccs", "repo"))
    touch(home / ".mzdesk/models/molgpka/molgpka_weights.npz")
    touch(home / ".mzdesk/models/neims/neims_weights.npz")
    assert reg["molgpka"].installed() and reg["neims"].installed()
    touch(home / ".mzdesk/data/ccs/ccs_index.json.gz")
    assert not reg["ccs"].installed()                                           # the meta file is missing too
    touch(home / ".mzdesk/data/ccs/ccs_meta.json")
    assert reg["ccs"].installed()
    touch(home / ".mzdesk/data/methods/report/report_index.npz")
    touch(home / ".mzdesk/data/methods/report/report_datasets.json")
    assert reg["repo"].installed()
    touch(home / ".mzdesk/data/libraries/MassBank/meta.json")
    assert not reg["libraries"].installed() and components.library_installed("MassBank") and not components.library_installed("GNPS")
    touch(home / ".mzdesk/data/libraries/GNPS/meta.json")
    assert reg["libraries"].installed()
    monkeypatch.setenv("MZDESK_CCS_DIR", str(home / "elsewhere"))
    assert not reg["ccs"].installed()                                           # the environment variable moves the folder


def test_msfinder_needs_mono_exe_databases_and_the_parameter_file(home, monkeypatch):
    folder = home / "msf"
    monkeypatch.setenv("MZDESK_MSFINDER_DIR", str(folder))
    monkeypatch.setenv("MZDESK_MONO", sys.executable)
    reg = components.REGISTRY
    assert not reg["msfinder"].installed()
    touch(folder / "MsfinderConsoleApp.exe")
    touch(folder / "Resources" / "MsfinderStructureDB-VS1.esd")
    assert not reg["msfinder"].installed()
    touch(folder / "Resources" / "MsfinderFormulaDB-VS1.efd")
    assert reg["msfinder"].installed()
    monkeypatch.setenv("MZDESK_MONO", str(home / "no-such-mono"))
    monkeypatch.setattr(paths, "mono_in_prefix", lambda: home / "none")
    monkeypatch.setattr("shutil.which", lambda *a, **k: None)
    from mzdesk.ms2 import msfinder_engine
    monkeypatch.setattr(msfinder_engine.os, "access", lambda *a: False)
    assert not reg["msfinder"].installed()                                      # no Mono, no MS-FINDER


def test_frigid_needs_every_part(home, monkeypatch):
    from mzdesk.ms2 import frigid_engine as fe
    folder = home / "frigid"
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(folder))
    comp = components.REGISTRY["frigid"]
    assert not comp.installed() and not fe.installed()
    py = touch(folder / "env/bin/python")
    py.chmod(0o755)
    touch(folder / "FRIGID/src/dlm/sampler.py", text="device = os.environ.get('FRIGID_DEVICE')")
    (folder / "FRIGID/src/mist").mkdir(parents=True)
    touch(folder / "env/lib/python3.10/site-packages/safe/__init__.py", text="# x\nfrom .converter import SAFEConverter, decode, encode\n")
    touch(folder / "checkpoints/DLM.ckpt", size=fe.DLM_SIZE)
    touch(folder / "checkpoints/mist_msg_aug.ckpt", size=fe.MIST_SIZE)
    assert not comp.installed()                                                 # the tokenizer is still missing
    touch(folder / "hf_cache/hub/models--datamol-io--safe-gpt/snapshots/abc/tokenizer.json")
    assert comp.installed() and fe.installed()
    (folder / "checkpoints/DLM.ckpt").write_text("half a download")
    assert not comp.installed()


def test_cfmid_follows_the_docker_state(monkeypatch):
    comp = components.REGISTRY["cfmid"]
    for state, expected in (("ready", True), ("no-image", False), ("down", False), ("no-docker", False)):
        monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True, s=state: s)
        assert comp.installed() is expected


def test_a_failing_status_check_means_not_installed():
    def broken():
        raise RuntimeError("boom")
    assert dataclasses.replace(components.REGISTRY["ccs"], check=broken).installed() is False


# ---------------------------------------------------------------- engine folders
def test_engine_folder_rule(home, monkeypatch):
    from mzdesk.ms2 import frigid_engine, msfinder_engine
    assert frigid_engine.install_dir() == home / ".mzdesk/tools/frigid"
    assert msfinder_engine.install_dir() == home / ".mzdesk/tools/msfinder/MSFINDER_v3.73"
    (home / "tools/frigid").mkdir(parents=True)
    (home / "tools/msfinder/MSFINDER_v3.73").mkdir(parents=True)
    assert frigid_engine.install_dir() == home / "tools/frigid"                                     # an existing ~/tools/... is kept
    assert msfinder_engine.install_dir() == home / "tools/msfinder/MSFINDER_v3.73"
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(home / "mine"))
    assert frigid_engine.install_dir() == home / "mine"                                             # the variable wins


def test_the_mono_made_by_the_installer_comes_before_the_one_on_path(home, monkeypatch):
    from mzdesk.ms2 import msfinder_engine
    own = touch(home / ".mzdesk/tools/mono/bin/mono")
    own.chmod(0o755)
    assert msfinder_engine.find_mono() == str(own)
    monkeypatch.setenv("MZDESK_MONO", sys.executable)
    assert msfinder_engine.find_mono() == sys.executable


def test_run_frigid_default_follows_the_install(home, monkeypatch):
    from mzdesk.ms2 import frigid_engine as fe
    monkeypatch.setattr(fe, "installed", lambda: False)
    assert fe.workflow_default() is False
    monkeypatch.setattr(fe, "installed", lambda: True)
    assert fe.workflow_default() is True
    monkeypatch.setenv("MZDESK_ANNOTATION_FRIGID", "off")
    assert fe.workflow_default() is False


def test_spec_reads_the_frigid_default_on_every_call(monkeypatch):
    from mzdesk import spec as spec_module
    from mzdesk.ms2 import frigid_engine as fe

    def run_frigid(body):
        pipe = next(p for p in body["pipelines"] if p["id"] == "esi_ms2_annotation")
        field = next(f for f in pipe["inputs"] if f["id"] == "run_frigid")
        return field, [e["values"]["run_frigid"] for e in pipe["examples"]]

    monkeypatch.setattr(fe, "installed", lambda: False)
    field, examples = run_frigid(spec_module.spec())
    assert field["needs"] == ["frigid"] and field["value"] is False and examples == [False, False]
    monkeypatch.setattr(fe, "installed", lambda: True)
    field, examples = run_frigid(spec_module.spec())
    assert field["value"] is True and examples == [True, True]
    assert [f for f in spec_module.PIPELINES[0]["inputs"] if f["id"] == "run_frigid"][0]["needs"] == ["frigid"]


# ---------------------------------------------------------------- needs()
def ids(items):
    return [i["id"] for i in items]


def test_every_module_and_workflow_target_is_known():
    from mzdesk.spec import MODULES, PIPELINES
    for m in MODULES:
        targets.needs(f"module:{m['id']}", {})
    for p in PIPELINES:
        targets.needs(f"workflow:{p['id']}", {})
    for bad in ("module:nope", "workflow:nope", "page:x", "", "ms2_esi"):
        with pytest.raises(targets.UnknownTarget):
            targets.needs(bad, {})


def test_module_targets(monkeypatch):
    with_installed(monkeypatch)
    expected = {"ms2_library": ["libraries"], "ms2_structure": ["msfinder"], "ms2_denovo": ["frigid"], "ms2_esi": [], "ms2_ei": ["neims"], "ccs": ["ccs"],
                "known_methods": ["repo", "molgpka"], "separation": ["molgpka"], "ionization": ["molgpka"], "adducts": ["molgpka"],
                "platform_fit": [], "non_ms": [], "derivatization": [], "ms2_fragments": [], "ms1_formula": [], "exclusivity": [], "identity": []}
    for module, want in expected.items():
        got = targets.needs(f"module:{module}", {})
        assert ids(got["missing"]) == want, module
        assert all(m["ask"] is components.REGISTRY[m["id"]].ask and m["title"] and m["size"] and m["why"] for m in got["missing"]), module
    with_installed(monkeypatch, msfinder=True, molgpka=True)
    assert targets.needs("module:ms2_structure", {}) == {"missing": [], "blocked": [], "unsupported": []}


def test_ms2_esi_is_blocked_by_docker_not_downloaded(monkeypatch):
    with_installed(monkeypatch)
    monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True: "no-docker")
    got = targets.needs("module:ms2_esi", {})
    assert got["missing"] == [] and got["blocked"][0]["id"] == "cfmid" and got["blocked"][0]["link"] == "https://docs.docker.com/get-docker/"
    monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True: "down")
    got = targets.needs("module:ms2_esi", {})
    assert got["blocked"][0]["link"] is None and "not running" in got["blocked"][0]["message"]
    monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True: "no-image")
    got = targets.needs("module:ms2_esi", {})
    assert got["blocked"] == [] and ids(got["missing"]) == ["cfmid"] and got["missing"][0]["ask"] is True
    monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True: "ready")
    assert targets.needs("module:ms2_esi", {})["missing"] == []


def test_library_search_asks_only_for_the_libraries_that_are_ticked(monkeypatch):
    with_installed(monkeypatch)
    monkeypatch.setattr(components, "library_installed", lambda name: name == "MassBank")
    assert ids(targets.needs("module:ms2_library", {"libraries": ["MassBank"]})["missing"]) == []
    assert ids(targets.needs("module:ms2_library", {"libraries": ["MassBank", "GNPS"]})["missing"]) == ["libraries"]
    assert ids(targets.needs("module:ms2_library", {"libraries": []})["missing"]) == []
    assert ids(targets.needs("module:ms2_library", {})["missing"]) == ["libraries"]            # no field: both libraries


def test_ms2_annotation_installs_libraries_and_msfinder_without_asking(monkeypatch):
    with_installed(monkeypatch)
    form = {"libraries": ["MassBank", "GNPS"], "run_frigid": False, "ccs_observed": ""}
    got = targets.needs("workflow:esi_ms2_annotation", form)
    assert [(m["id"], m["ask"]) for m in got["missing"]] == [("libraries", False), ("msfinder", False), ("molgpka", False)]
    # the same components on their own pages do ask
    assert targets.needs("module:ms2_library", {})["missing"][0]["ask"] is True
    assert targets.needs("module:ms2_structure", {})["missing"][0]["ask"] is True


def test_ms2_annotation_frigid_and_ccs_conditions(monkeypatch):
    with_installed(monkeypatch)
    base = {"libraries": [], "run_frigid": False}
    assert ids(targets.needs("workflow:esi_ms2_annotation", base)["missing"]) == ["msfinder", "molgpka"]                 # no libraries ticked: none needed
    got = targets.needs("workflow:esi_ms2_annotation", {**base, "run_frigid": True})
    assert [(m["id"], m["ask"]) for m in got["missing"]][-1] == ("frigid", True)                                       # FRIGID always asks
    assert "frigid" not in ids(targets.needs("workflow:esi_ms2_annotation", {**base, "run_frigid": "false"})["missing"])
    for observed, wanted in (("", False), (None, False), ("  ", False), ("152.3", True), (152.3, True)):
        got = targets.needs("workflow:esi_ms2_annotation", {**base, "ccs_observed": observed})
        assert ("ccs" in ids(got["missing"])) is wanted, observed
    ccs = next(m for m in targets.needs("workflow:esi_ms2_annotation", {**base, "ccs_observed": "152.3"})["missing"] if m["id"] == "ccs")
    assert ccs["ask"] is False
    with_installed(monkeypatch, msfinder=True, molgpka=True, ccs=True, frigid=True)
    assert targets.needs("workflow:esi_ms2_annotation", {**base, "run_frigid": True, "ccs_observed": "152.3"})["missing"] == []


def test_other_workflows(monkeypatch):
    with_installed(monkeypatch)
    got = targets.needs("workflow:ei_ms_annotation", {})
    assert [(m["id"], m["ask"]) for m in got["missing"]] == [("libraries", True), ("neims", True)]
    got = targets.needs("workflow:targeted_method", {})
    assert ids(got["missing"]) == ["molgpka", "repo"] and not any(m["ask"] for m in got["missing"])


def test_repo_brings_molgpka_along(monkeypatch):
    with_installed(monkeypatch)
    monkeypatch.setattr(targets, "MODULES", {**targets.MODULES, "known_methods": ("repo",)})
    assert ids(targets.needs("module:known_methods", {})["missing"]) == ["repo", "molgpka"]


def test_unsupported_components_are_listed_and_do_not_block(monkeypatch):
    with_installed(monkeypatch)
    monkeypatch.setattr(components, "current_platform", lambda: "win32")
    got = targets.needs("workflow:esi_ms2_annotation", {"libraries": ["GNPS"], "run_frigid": True, "ccs_observed": "150"})
    assert ids(got["missing"]) == ["libraries", "molgpka", "ccs"] and ids(got["unsupported"]) == ["msfinder", "frigid"]
    assert got["unsupported"][0]["message"] == "macOS / Linux only for now" and got["blocked"] == []
    assert targets.needs("module:ms2_esi", {})["unsupported"][0]["id"] == "cfmid"


# ---------------------------------------------------------------- install jobs
def fake_component(cid, run, state, **kw):
    def check():
        return state.get(cid, False)

    def runner(ctx):
        run(ctx)
        if not kw.get("leave_incomplete"):
            state[cid] = True
    return components.Component(cid, f"Fake {cid}", "1 MB", kw.get("ask", False), components.ALL_PLATFORMS, check, runner, refresh=kw.get("refresh"))


def child_script(code):
    return [sys.executable, "-u", "-c", code]


def wait(job, seconds=30):
    job.thread.join(seconds)
    assert job.finished(), job.snapshot()
    return job.snapshot()


@pytest.fixture
def fakes(monkeypatch):
    state = {}
    steps = []

    def one(ctx):
        steps.append("one")
        ctx.run(child_script("print('Downloading a'); print('50 %', end='\\r'); print('done a')"), what="step a")

    def two(ctx):
        steps.append("two")
        ctx.log("second component")

    reg = {"one": fake_component("one", one, state), "two": fake_component("two", two, state)}
    monkeypatch.setattr(components, "REGISTRY", reg)
    return state, steps


def test_job_runs_components_in_order_and_streams_the_log(fakes):
    state, steps = fakes
    store = jobs.JobStore()
    job = store.start(["two", "one"])
    snap = wait(job)
    assert steps == ["one", "two"] and snap["status"] == "done" and snap["error"] == "" and snap["current"] is None
    assert [(c["id"], c["status"]) for c in snap["components"]] == [("one", "done"), ("two", "done")]
    assert snap["log"][:3] == ["Downloading a", "50 %", "done a"] and snap["log"][-1] == "second component"
    assert snap["elapsed_s"] >= 0 and set(snap) == {"status", "current", "elapsed_s", "components", "log", "error"}
    assert set(snap["components"][0]) == {"id", "title", "status", "detail"} and snap["components"][1]["detail"] == "installed"
    assert store.get(job.id) is job and store.get("nope") is None


def test_job_keeps_only_the_last_30_lines(monkeypatch):
    state = {}
    reg = {"big": fake_component("big", lambda ctx: [ctx.log(f"line {i}") for i in range(100)], state)}
    monkeypatch.setattr(components, "REGISTRY", reg)
    snap = wait(jobs.JobStore().start(["big"]))
    assert len(snap["log"]) == 30 and snap["log"][-1] == "line 99"


def test_a_component_that_is_installed_is_skipped(fakes):
    state, steps = fakes
    state["one"] = True
    snap = wait(jobs.JobStore().start(["one", "two"]))
    assert steps == ["two"] and snap["components"][0]["detail"] == "already installed" and snap["status"] == "done"


def test_a_failing_step_ends_the_job_with_the_message_and_stops_the_rest(monkeypatch):
    state, steps = {}, []

    def broken(ctx):
        ctx.run(child_script("print('about to fail'); raise SystemExit(3)"), what="the broken step")

    reg = {"a": fake_component("a", broken, state), "b": fake_component("b", lambda ctx: steps.append("b"), state)}
    monkeypatch.setattr(components, "REGISTRY", reg)
    snap = wait(jobs.JobStore().start(["a", "b"]))
    assert snap["status"] == "error" and "the broken step failed (exit code 3)" in snap["error"] and "about to fail" in snap["error"]
    assert [c["status"] for c in snap["components"]] == ["error", "waiting"] and steps == []
    assert snap["components"][0]["detail"].startswith("the broken step failed")


def test_unexpected_exception_and_incomplete_install_are_errors(monkeypatch):
    state = {}

    def boom(ctx):
        raise KeyError("x")

    reg = {"a": fake_component("a", boom, state), "b": fake_component("b", lambda ctx: None, state, leave_incomplete=True)}
    monkeypatch.setattr(components, "REGISTRY", reg)
    snap = wait(jobs.JobStore().start(["a"]))
    assert snap["status"] == "error" and "KeyError" in snap["error"]
    snap = wait(jobs.JobStore().start(["b"]))
    assert snap["status"] == "error" and "not complete" in snap["error"]


def test_cancel_kills_the_child_process_tree(monkeypatch, home):
    state = {}
    marker = home / "grandchild.pid"
    code = ("import subprocess, sys, time\n"
            f"p = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(120)'])\n"
            f"open({str(marker)!r}, 'w').write(str(p.pid))\n"
            "print('started', flush=True)\n"
            "time.sleep(120)\n")
    reg = {"slow": fake_component("slow", lambda ctx: ctx.run(child_script(code), what="slow"), state), "after": fake_component("after", lambda ctx: None, state)}
    monkeypatch.setattr(components, "REGISTRY", reg)
    store = jobs.JobStore()
    job = store.start(["slow", "after"])
    deadline = time.monotonic() + 20
    while not marker.exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    assert marker.exists() and job.snapshot()["current"] == "slow"
    assert store.cancel(job.id) is True
    snap = wait(job, 30)
    assert snap["status"] == "cancelled" and snap["error"] == "" and not state
    pid = int(marker.read_text())
    time.sleep(0.3)
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)                                         # the grandchild is gone too
    assert store.cancel("unknown") is False


def test_only_one_job_runs_at_a_time(monkeypatch):
    state = {}
    gate = threading.Event()
    reg = {"hold": fake_component("hold", lambda ctx: gate.wait(20), state), "other": fake_component("other", lambda ctx: None, state)}
    monkeypatch.setattr(components, "REGISTRY", reg)
    store = jobs.JobStore()
    first = store.start(["hold"])
    assert store.start(["other"]) is first                       # a second call returns the running job
    gate.set()
    wait(first)
    assert store.start(["other"]) is not first


def test_old_jobs_are_pruned(monkeypatch):
    state = {}
    monkeypatch.setattr(components, "REGISTRY", {"a": fake_component("a", lambda ctx: None, state)})
    store = jobs.JobStore()
    first = store.start(["a"])
    wait(first)
    state.clear()
    for _ in range(jobs.KEEP_JOBS + 3):
        wait(store.start(["a"]))
        state.clear()
    assert len(store.jobs) <= jobs.KEEP_JOBS and first.id not in store.jobs
    monkeypatch.setattr(jobs, "TTL_S", -1)
    store.start(["a"])
    assert len(store.jobs) <= 2


def test_refresh_runs_after_a_successful_install(monkeypatch):
    state, called = {}, []
    monkeypatch.setattr(components, "REGISTRY", {"a": fake_component("a", lambda ctx: None, state, refresh=lambda: called.append(1))})
    wait(jobs.JobStore().start(["a"]))
    assert called == [1]


def test_unknown_or_unsupported_ids_do_not_start_a_job(fakes):
    with pytest.raises(InstallError):
        jobs.JobStore().start(["nope"])


# ---------------------------------------------------------------- ctx
def test_ctx_splits_on_carriage_returns_strips_ansi_and_reports_the_exit_code():
    lines = []
    ctx = Ctx(lines.append)
    code = ctx.run(child_script("import sys; sys.stdout.write('a\\rb\\n\\x1b[2Kc\\n'); sys.stdout.flush(); sys.stderr.write('err\\n')"), what="x", check=False)
    assert code == 0 and lines[:2] == ["a", "b"] and "c" in lines and "err" in lines
    with pytest.raises(InstallError, match="x failed .exit code 2"):
        ctx.run(child_script("raise SystemExit(2)"), what="x")
    with pytest.raises(InstallError, match="Could not start"):
        ctx.run(["/no/such/program"])
    ctx.cancel()
    with pytest.raises(Cancelled):
        ctx.run(child_script("print(1)"))


# ---------------------------------------------------------------- commit archive without git
def make_archive(commit, files):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for name, data in files.items():
            info = tarfile.TarInfo(f"repo-{commit}/{name}")
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
        link = tarfile.TarInfo(f"repo-{commit}/link")
        link.type, link.linkname = tarfile.SYMTYPE, "/etc/passwd"
        tf.addfile(link)
    return buf.getvalue()


@pytest.fixture
def archive_server(monkeypatch):
    served = {}

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            body = served.get("body")
            self.send_response(200 if body else 404)
            self.end_headers()
            self.wfile.write(body or b"")

        def log_message(self, *a):
            pass

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    monkeypatch.setattr(archive, "archive_url", lambda repo_url, commit: f"http://127.0.0.1:{server.server_port}/{commit}.tar.gz")
    yield served
    server.shutdown()


COMMIT = "a" * 40
PATTERNS = ["/LICENSE.txt", "/data/*.json", "/processed/*/*_info.tsv"]


def test_archive_keeps_only_the_wanted_files(archive_server, tmp_path):
    archive_server["body"] = make_archive(COMMIT, {"LICENSE.txt": b"mit", "README.md": b"r", "data/a.json": b"{}", "data/sub/b.json": b"{}", "data/c.py": b"x",
                                                   "processed/1/1_info.tsv": b"t", "processed/1/1_other.tsv": b"t", "sub/LICENSE.txt": b"no"})
    dest = tmp_path / "out"
    n = archive.fetch("https://github.com/o/r", COMMIT, dest, PATTERNS, log=lambda text: None)
    files = sorted(p.relative_to(dest).as_posix() for p in dest.rglob("*") if p.is_file())
    assert n == 3 and files == [archive.MARKER, "LICENSE.txt", "data/a.json", "processed/1/1_info.tsv"]       # no symlink, no deeper '*', no sibling names
    assert archive.has_commit(dest, COMMIT) and not archive.has_commit(dest, "b" * 40) and not archive.has_commit(tmp_path / "none", COMMIT)


def test_archive_errors_are_clear(archive_server, tmp_path):
    with pytest.raises(archive.ArchiveError, match="Download failed"):
        archive.fetch("https://github.com/o/r", COMMIT, tmp_path / "a", PATTERNS, log=print)               # 404
    archive_server["body"] = make_archive("b" * 40, {"LICENSE.txt": b"x"})
    with pytest.raises(archive.ArchiveError, match="expected a folder that ends with"):
        archive.fetch("https://github.com/o/r", COMMIT, tmp_path / "b", PATTERNS, log=print)
    archive_server["body"] = make_archive(COMMIT, {"nothing.txt": b"x"})
    with pytest.raises(archive.ArchiveError, match="No file of the archive matched"):
        archive.fetch("https://github.com/o/r", COMMIT, tmp_path / "c", PATTERNS, log=print)
    archive_server["body"] = b"not a tar file at all"
    with pytest.raises(archive.ArchiveError, match="interrupted or damaged"):
        archive.fetch("https://github.com/o/r", COMMIT, tmp_path / "d", PATTERNS, log=print)


def test_archive_url_and_patterns():
    assert archive.archive_url("https://github.com/dylanhross/c3sdb", COMMIT) == f"https://codeload.github.com/dylanhross/c3sdb/tar.gz/{COMMIT}"
    assert archive.archive_url("https://github.com/michaelwitting/RepoRT.git", COMMIT).endswith(f"/michaelwitting/RepoRT/tar.gz/{COMMIT}")
    for bad in (("https://evil.example/o/r", COMMIT), ("https://github.com/o/r", "main"), ("https://github.com/o", COMMIT)):
        with pytest.raises(archive.ArchiveError):
            archive.archive_url(*bad)
    assert archive.pattern_regex("/processed_data/*/*_info.tsv").fullmatch("processed_data/0001/0001_info.tsv")
    assert not archive.pattern_regex("/processed_data/*/*_info.tsv").fullmatch("processed_data/a/b/0001_info.tsv")


def test_ccs_and_repo_setups_use_the_archive_without_git(archive_server, tmp_path, monkeypatch):
    from mzdesk.methods import report_setup
    from mzdesk.orthogonal import ccs_setup
    monkeypatch.setattr(paths, "working_git", lambda: None)
    archive_server["body"] = make_archive(ccs_setup.COMMIT, {"LICENSE.txt": b"MIT", "README.md": b"r", "c3sdb/_include/smiles_search_cache.json": b"{}",
                                                             "c3sdb/_include/src_data/x.json": b"{}", "c3sdb/other.py": b"no"})
    ccs_setup.fetch(tmp_path / "ccs", ccs_setup.COMMIT)
    got = sorted(p.relative_to(tmp_path / "ccs").as_posix() for p in (tmp_path / "ccs").rglob("*") if p.is_file())
    assert got == [archive.MARKER, "LICENSE.txt", "README.md", "c3sdb/_include/smiles_search_cache.json", "c3sdb/_include/src_data/x.json"]
    ccs_setup.fetch(tmp_path / "ccs", ccs_setup.COMMIT)                                                    # second call: nothing to download
    archive_server["body"] = make_archive(report_setup.COMMIT, {"LICENSE": b"x", "processed_data/0001/0001_info.tsv": b"t", "processed_data/0001/0001_x.bin": b"b"})
    report_setup.fetch(tmp_path / "report", report_setup.COMMIT)
    got = sorted(p.relative_to(tmp_path / "report").as_posix() for p in (tmp_path / "report").rglob("*") if p.is_file())
    assert got == [archive.MARKER, "LICENSE", "processed_data/0001/0001_info.tsv"]


# ---------------------------------------------------------------- Java and git discovery
def test_java_on_path_must_start_else_jdk4py(monkeypatch, tmp_path):
    paths.find_java.cache_clear()
    fake_java = tmp_path / "java"
    fake_java.write_text("#!/bin/sh\necho 'Unable to locate a Java Runtime' >&2\nexit 1\n")
    fake_java.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda name, *a, **k: str(fake_java) if name == "java" else None)
    monkeypatch.setitem(sys.modules, "jdk4py", None)                                                       # import fails
    assert paths.find_java() is None
    paths.find_java.cache_clear()
    jdk = types_module("jdk4py", JAVA=tmp_path / "bundled-java")
    (tmp_path / "bundled-java").write_text("")
    monkeypatch.setitem(sys.modules, "jdk4py", jdk)
    assert paths.find_java() == str(tmp_path / "bundled-java")                                             # the stub on PATH is refused
    paths.find_java.cache_clear()
    fake_java.write_text("#!/bin/sh\nexit 0\n")
    assert paths.find_java() == str(fake_java)
    paths.find_java.cache_clear()


def types_module(name, **attrs):
    import types
    mod = types.ModuleType(name)
    mod.__dict__.update(attrs)
    return mod


def test_metfrag_reports_the_missing_java(monkeypatch):
    from mzdesk.ms2 import metfrag_engine as mf
    monkeypatch.setattr(paths, "find_java", lambda: None)
    with pytest.raises(mf.MetFragError, match="jdk4py"):
        mf.run_metfrag("CC(C)O", [(43.0542, 100.0)], "[M+H]+", None, 2, 0.01, 0.0)


def test_working_git_ignores_the_macos_stub_without_developer_tools(monkeypatch, tmp_path):
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.setattr("shutil.which", lambda name, *a, **k: "/usr/bin/git")
    monkeypatch.setattr(Path, "is_file", lambda self: str(self) == "/usr/bin/git")

    class Done:
        returncode = 2
    monkeypatch.setattr(paths.subprocess, "run", lambda *a, **k: Done())
    assert paths.working_git() is None
    Done.returncode = 0
    assert paths.working_git() == "/usr/bin/git"


# ---------------------------------------------------------------- prerequisites that need no network
def test_micromamba_platform_names(monkeypatch):
    for system, machine, expected in (("darwin", "arm64", "osx-arm64"), ("darwin", "x86_64", "osx-64"), ("linux", "x86_64", "linux-64"),
                                      ("linux", "aarch64", "linux-aarch64"), ("win32", "AMD64", "win-64")):
        monkeypatch.setattr(sys, "platform", system)
        monkeypatch.setattr(prereq.platform, "machine", lambda m=machine: m)
        assert prereq.micromamba_platform() == expected
    monkeypatch.setattr(prereq.platform, "machine", lambda: "riscv64")
    with pytest.raises(InstallError, match="no micromamba build"):
        prereq.micromamba_platform()


def test_a_mono_that_works_is_used_without_downloading(home, monkeypatch):
    own = touch(home / ".mzdesk/tools/mono/bin/mono")
    own.chmod(0o755)
    monkeypatch.setattr(prereq, "mono_works", lambda mono: (True, ""))
    monkeypatch.setattr(prereq, "conda_create", lambda *a, **k: pytest.fail("must not create an environment"))
    mono, own_bin = prereq.ensure_mono(Ctx(lambda line: None))
    assert mono == str(own) and own_bin == own.parent


def test_mono_on_apple_silicon_needs_rosetta_and_an_osx_64_environment(home, monkeypatch):
    monkeypatch.setattr(prereq, "mono_works", lambda mono: (False, "no csc"))
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.setattr(prereq.platform, "machine", lambda: "arm64")
    monkeypatch.setattr(prereq, "rosetta_ready", lambda: False)
    with pytest.raises(InstallError, match="softwareupdate --install-rosetta"):
        prereq.ensure_mono(Ctx(lambda line: None), system=False)
    made = {}

    def fake_create(ctx, prefix, specs, subdir=None, what=""):
        made.update(prefix=prefix, specs=specs, subdir=subdir)
        mono = Path(prefix) / "bin" / "mono"
        mono.parent.mkdir(parents=True)
        mono.write_text("")
        mono.chmod(0o755)
    monkeypatch.setattr(prereq, "rosetta_ready", lambda: True)
    monkeypatch.setattr(prereq, "conda_create", fake_create)
    monkeypatch.setattr(prereq, "mono_works", lambda mono: (Path(mono).is_file() and "tools/mono" in mono, "x"))
    mono, own_bin = prereq.ensure_mono(Ctx(lambda line: None), system=False)
    assert made["subdir"] == "osx-64" and made["specs"] == ["mono=6.12.0.199"] and made["prefix"] == home / ".mzdesk/tools/mono" and own_bin == Path(mono).parent


def test_throwaway_environment_command_never_uses_the_main_environment(home, monkeypatch):
    seen = {}
    monkeypatch.setattr(prereq, "find_uv", lambda: "/fake/uv")
    monkeypatch.setattr(sys, "platform", "linux")

    class FakeCtx:
        def log(self, text):
            pass

        def run(self, cmd, env=None, **kw):
            seen.update(cmd=cmd, env=env)
            assert Path(env["UV_CACHE_DIR"]).is_dir()

    monkeypatch.setenv("VIRTUAL_ENV", "/some/venv")
    monkeypatch.setenv("PYTHONPATH", "/leak")
    prereq.throwaway_python(FakeCtx(), ["torch", "numpy"], "/x/pka_setup.py", "the conversion")
    cmd, env = seen["cmd"], seen["env"]
    assert cmd[:2] == ["/fake/uv", "run"] and "--no-project" in cmd and cmd[-3:] == ["python", "-I", "/x/pka_setup.py"]
    assert [cmd[i + 1] for i, c in enumerate(cmd) if c == "--with"] == ["torch", "numpy"] and "https://download.pytorch.org/whl/cpu" in cmd
    assert "VIRTUAL_ENV" not in env and "PYTHONPATH" not in env and not Path(env["UV_CACHE_DIR"]).exists()      # the cache is deleted afterwards


def test_missing_uv_says_how_to_get_it(monkeypatch):
    monkeypatch.setattr(prereq, "find_uv", lambda: None)
    with pytest.raises(InstallError, match="uv was not found.*install.sh"):
        prereq.throwaway_python(Ctx(lambda line: None), ["torch"], "/x.py", "x")


def test_uv_is_looked_up_in_mzdesk_bin_then_on_path(home, monkeypatch):
    on_path = touch(home / "elsewhere/uv")
    on_path.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda name, *a, **k: str(on_path) if name == "uv" else None)
    assert prereq.find_uv() == str(on_path)                                      # the installer may have left uv on PATH only
    own = touch(home / ".mzdesk/bin/uv")
    own.chmod(0o755)
    assert prereq.find_uv() == str(own)                                          # ~/.mzdesk/bin/uv comes first
    own.unlink()
    on_path.unlink()
    local = touch(home / ".local/bin/uv")
    local.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda *a, **k: None)
    assert prereq.find_uv() == str(local)
    local.unlink()
    assert prereq.find_uv() is None
