"""Windows paths of the local engines. Most tests run anywhere (they fake the platform or the files); the ones marked `windows` need real Windows
and run in the Windows CI job (.github/workflows/windows.yml)."""
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from mzdesk.ms2 import cfmid_engine as cf
from mzdesk.ms2 import frigid_engine as fe
from mzdesk.ms2 import frigid_setup as fs
from mzdesk.ms2 import msfinder_engine as mf
from mzdesk.ms2 import msfinder_setup as ms
from mzdesk.shared import memguard, procgroup

windows = pytest.mark.skipif(sys.platform != "win32", reason="needs Windows")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------- MS-FINDER without Mono
def test_the_bundled_windows_files_are_the_pinned_ones():
    assert sha(ms.WIN_EXE) == ms.WIN_EXE_SHA256
    assert sha(ms.WIN_FIXED_DLL) == mf.STRUCTUREFINDER_FIXED_SHA256
    assert (ms.WIN_DIR / "NOTICE.md").is_file()


def test_launch_command_runs_the_exe_alone_on_windows_and_through_mono_elsewhere():
    assert mf.launch_command(mf.NATIVE, Path("C:/x/MsfinderConsoleApp.exe")) == [str(Path("C:/x/MsfinderConsoleApp.exe"))]
    assert mf.launch_command("/usr/bin/mono", "/x/a.exe") == ["/usr/bin/mono", "/x/a.exe"]


def test_find_mono_on_windows_is_native_when_the_framework_is_recent(monkeypatch, tmp_path):
    monkeypatch.setattr(mf.sys, "platform", "win32")
    monkeypatch.setenv("WINDIR", str(tmp_path))
    assert mf.find_mono() is None
    framework = tmp_path / "Microsoft.NET" / "Framework64" / "v4.0.30319"
    framework.mkdir(parents=True)
    (framework / "netstandard.dll").write_bytes(b"")
    assert mf.find_mono() == mf.NATIVE


def test_install_prebuilt_makes_the_same_end_state_as_the_mono_build(monkeypatch, tmp_path):
    release = b"pretend release StructureFinder.dll"
    install = tmp_path / "MSFINDER_v3.73"
    install.mkdir()
    (install / mf.STRUCTUREFINDER_DLL).write_bytes(release)
    monkeypatch.setattr(mf, "STRUCTUREFINDER_SHA256", hashlib.sha256(release).hexdigest())
    assert mf.thread_fix_state(install) == "release"
    assert ms.install_prebuilt(install) == "fixed"
    assert sha(install / mf.EXE_NAME) == ms.WIN_EXE_SHA256
    assert (install / (mf.EXE_NAME + ".config")).read_bytes() == mf.WIN_CONFIG.read_bytes()
    assert (install / (mf.STRUCTUREFINDER_DLL + ".original")).read_bytes() == release
    marker = json.loads((install / mf.THREAD_FIX_MARKER).read_text(encoding="utf-8"))
    assert marker["fixed_sha256"] == mf.STRUCTUREFINDER_FIXED_SHA256 and marker["same_as_reference_build"] is True
    assert ms.install_prebuilt(install) == "fixed"                          # again: nothing changes, nothing fails


def test_the_console_exe_config_redirects_to_the_release_builds():
    import xml.etree.ElementTree as ET
    ns = {"a": "urn:schemas-microsoft-com:asm.v1"}
    root = ET.parse(mf.WIN_CONFIG).getroot()
    got = {d.find("a:assemblyIdentity", ns).get("name"): (d.find("a:assemblyIdentity", ns).get("publicKeyToken"), d.find("a:bindingRedirect", ns).get("newVersion"))
           for d in root.iter("{urn:schemas-microsoft-com:asm.v1}dependentAssembly")}
    # versions and tokens of the DLLs in the MS-FINDER v3.73 release folder; MessagePack asks for System.Threading.Tasks.Extensions 4.1.0.0
    assert got == {"System.Threading.Tasks.Extensions": ("cc7b13ffcd2ddd51", "4.2.0.1"), "System.Buffers": ("cc7b13ffcd2ddd51", "4.0.3.0"),
                   "System.Memory": ("cc7b13ffcd2ddd51", "4.0.1.2"), "System.Runtime.CompilerServices.Unsafe": ("b03f5f7f11d50a3a", "4.0.4.1"),
                   "System.Numerics.Vectors": ("b03f5f7f11d50a3a", "4.1.4.0"), "System.ValueTuple": ("cc7b13ffcd2ddd51", "4.0.3.0")}


@pytest.mark.parametrize("platform", ["win32", "darwin"])
def test_each_run_copy_has_the_config_on_windows_only(monkeypatch, tmp_path, platform):
    src = tmp_path / "install"
    src.mkdir()
    (src / mf.EXE_NAME).write_bytes(b"exe")
    monkeypatch.setattr(mf.sys, "platform", platform)
    monkeypatch.setattr(mf, "_clone_tree", lambda a, b: (_ for _ in ()).throw(OSError("no clone here")))  # plain copy: the real clone needs the real platform
    exe, _ = mf.private_install(src / mf.EXE_NAME, tmp_path / "run")
    config = exe.parent / (mf.EXE_NAME + ".config")
    assert config.is_file() == (platform == "win32")
    assert not (src / (mf.EXE_NAME + ".config")).exists()                   # the shared install is only read
    if platform == "win32":
        assert config.read_bytes() == mf.WIN_CONFIG.read_bytes()


def test_install_prebuilt_refuses_a_dll_that_is_not_the_release(tmp_path):
    install = tmp_path / "MSFINDER_v3.73"
    install.mkdir()
    (install / mf.STRUCTUREFINDER_DLL).write_bytes(b"something else")
    with pytest.raises(ms.SetupError, match="not the MS-FINDER v3.73 release"):
        ms.install_prebuilt(install)


def test_install_prebuilt_refuses_a_damaged_bundle(monkeypatch, tmp_path):
    bad = tmp_path / "bad.exe"
    bad.write_bytes(b"x")
    monkeypatch.setattr(ms, "WIN_EXE", bad)
    with pytest.raises(ms.SetupError, match="expected"):
        ms.install_prebuilt(tmp_path)


# ---------------------------------------------------------------- FRIGID
def test_env_python_and_safe_init_follow_the_conda_layout_of_windows(monkeypatch, tmp_path):
    monkeypatch.setattr(fe.sys, "platform", "linux")
    assert fe.env_python(tmp_path) == tmp_path / "env" / "bin" / "python"
    monkeypatch.setattr(fe.sys, "platform", "win32")
    assert fe.env_python(tmp_path) == tmp_path / "env" / "python.exe"
    init = tmp_path / "env" / "Lib" / "site-packages" / "safe" / "__init__.py"
    init.parent.mkdir(parents=True)
    init.write_text("x\n", encoding="utf-8")
    assert fe.safe_init(tmp_path) == init


def test_conda_candidates_on_windows_use_the_scripts_folder(monkeypatch):
    monkeypatch.setattr(fs.sys, "platform", "win32")
    assert all(c.endswith("conda.exe") and "Scripts" in c for c in fs.conda_candidates())


def test_safe_patch_writes_unix_line_ends(tmp_path):
    env = tmp_path / "env" / "lib" / "python3.10" / "site-packages" / "safe"
    env.mkdir(parents=True)
    (env / "__init__.py").write_text("import a\nimport b\n", encoding="utf-8", newline="\n")
    assert fs.apply_safe_patch(tmp_path) == "applied"
    raw = (env / "__init__.py").read_bytes()
    assert b"\r" not in raw and raw.endswith(fs.SAFE_INIT_LAST_LINE.encode() + b"\n")
    assert fs.apply_safe_patch(tmp_path) == "already applied"


# ---------------------------------------------------------------- CFM-ID
def test_cfmid_mounts_the_work_folder_with_mount_not_v(monkeypatch, tmp_path):
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")
    monkeypatch.setattr(cf.subprocess, "run", fake_run)
    cf._docker_run(["echo"], str(tmp_path))
    cmd = seen["cmd"]
    assert "-v" not in cmd and f"type=bind,source={tmp_path},target=/data" in cmd


# ---------------------------------------------------------------- process groups (runs on every system; the Windows branch in CI)
def _sleeper():
    return subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"], **procgroup.popen_kwargs())


def test_a_child_in_its_own_group_can_be_measured_and_killed():
    proc = _sleeper()
    try:
        time.sleep(0.5)
        assert procgroup.group_rss_mb(proc.pid) > 1
        procgroup.kill_group(proc)
        assert proc.poll() is not None
    finally:
        if proc.poll() is None:
            proc.kill()


@windows
def test_memguard_reads_windows_memory():
    assert memguard.system_memory()["total_mb"] > 1000
    mine = memguard.proc_memory(os.getpid())
    assert mine and mine["rss_mb"] > 10
    proc = _sleeper()
    try:
        time.sleep(0.5)
        assert proc.pid in memguard.tree_pids(os.getpid())
        assert memguard.tree_memory(os.getpid())["n_procs"] >= 2
    finally:
        proc.kill()


def test_pid_alive_tells_a_running_child_from_a_gone_one():
    proc = _sleeper()
    try:
        assert procgroup.pid_alive(proc.pid)
        procgroup.kill_group(proc)
        proc.wait(timeout=10)
        time.sleep(0.2)
        assert not procgroup.pid_alive(proc.pid)
    finally:
        if proc.poll() is None:
            proc.kill()


@windows
def test_msfinder_framework_is_there_on_a_current_windows():
    assert mf.find_mono() == mf.NATIVE
