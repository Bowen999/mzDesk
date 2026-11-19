"""GET /api/versions and mzdesk/shared/versions.py: the versions the Session report writes into report.json.

Nothing here starts Docker, Mono or Java: the three fixed commands go through versions._run, which the tests replace.
"""
import json
import os
import platform
import re
import sys
import types
from datetime import datetime
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk import __version__
from mzdesk.api.app import app
from mzdesk.ms2 import cfmid_engine
from mzdesk.shared import versions
from mzdesk.spec import MODULES, PIPELINES

client = TestClient(app)
ENTRY_KEYS = {"version", "kind", "detail"}
ABSOLUTE_PATH = re.compile(r"(?<![\w.:])/(?:Users|home|opt|var|usr|tmp|etc|private|srv|mnt|root|Library|Applications)/|[A-Za-z]:\\")


@pytest.fixture(autouse=True)
def fresh():
    """Every test starts without a cached answer and leaves none behind."""
    versions.reset()
    for probe in (versions._docker_image, versions._mono_version, versions._java_version):
        probe.cache_clear()
    yield
    versions.reset()


@pytest.fixture
def nothing_installed(tmp_path, monkeypatch):
    """Every folder the engines read points at an empty place, and no command can run."""
    empty = tmp_path / "empty"
    empty.mkdir()
    for var in ("MZDESK_MSFINDER_DIR", "MZDESK_FRIGID_DIR", "MZDESK_LIBRARY_DIR", "MZDESK_CCS_DIR", "MZDESK_METHODS_DIR"):
        monkeypatch.setenv(var, str(empty / var.lower()))
    for var in ("MZDESK_REPORT_DIR", "MZDESK_MONO"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("MZDESK_NEIMS_WEIGHTS", str(empty / "neims.npz"))
    monkeypatch.setenv("MZDESK_MOLGPKA_WEIGHTS", str(empty / "molgpka.npz"))
    monkeypatch.setattr(versions.shutil, "which", lambda name, *a, **k: None)
    monkeypatch.setattr(versions, "_run", lambda cmd, timeout: None)
    return empty


def based_on_names():
    """Every name a page lists under 'Based on' (all pages, also the ones held back)."""
    names = set()
    for page in [*MODULES, *PIPELINES]:
        names.update(page.get("based_on") or [])
    return names


# ---------------------------------------------------------------- the endpoint
def test_versions_endpoint_has_the_documented_shape():
    r = client.get("/api/versions")
    assert r.status_code == 200
    body = r.json()
    assert set(body) == {"mzdesk", "python", "checked_at", "components"}
    assert body["mzdesk"] == __version__
    assert body["python"] == platform.python_version()
    assert datetime.fromisoformat(body["checked_at"]).tzinfo is not None
    assert body["components"]["mzDesk"] == {"version": __version__, "kind": "mzdesk", "detail": None}
    for name, entry in body["components"].items():
        assert set(entry) == ENTRY_KEYS, name
        assert entry["kind"] in versions.KINDS, name
        assert entry["version"] is None or (isinstance(entry["version"], str) and entry["version"]), name
        assert entry["detail"] is None or (isinstance(entry["detail"], str) and entry["detail"]), name
        if entry["version"] is None:
            assert entry["detail"], f"{name}: a missing version says why"
    json.dumps(body)


def test_spec_is_unchanged_by_the_versions_endpoint():
    assert {"groups", "pipelines", "subgroups", "modules", "version", "hosted"} <= set(client.get("/api/spec").json())


def test_every_based_on_name_of_every_page_has_a_resolver():
    names = based_on_names()
    assert names, "the pages list what they are based on"
    assert names <= set(versions.RESOLVERS), sorted(names - set(versions.RESOLVERS))
    assert names <= set(client.get("/api/versions").json()["components"])
    assert "GNPS (online)" in versions.RESOLVERS, "build.ts reads it for the pages that use GNPS online"


def test_online_services_are_marked_online_and_have_no_version():
    components = client.get("/api/versions").json()["components"]
    for name in ("fastMASST", "Pan-ReDU", "PubChem", "UniChem", "LOTUS", "ClassyFire", "GNPS (online)"):
        assert components[name]["kind"] == "online" and components[name]["version"] is None, name
        assert "not versioned" in components[name]["detail"]


def test_key_packages_resolve_to_the_installed_versions():
    components = client.get("/api/versions").json()["components"]
    from importlib.metadata import distributions
    for name, dist in (("RDKit", "rdkit"), ("numpy", "numpy"), ("BUDDY", "msbuddy")):      # one of the versions the metadata lists (an upgrade can leave two)
        assert components[name]["version"] in {d.version for d in distributions(name=dist)}, name
        assert components[name]["kind"] == "python-package"


# ---------------------------------------------------------------- a missing thing is null, never an error
def test_a_missing_package_gives_a_null_version_and_does_not_raise():
    assert versions.python_package("no-such-package-xyz") == {"version": None, "kind": "python-package", "detail": "not installed"}


def test_broken_metadata_does_not_raise(monkeypatch):
    def boom(**kwargs):
        raise RuntimeError("broken metadata folder")
    monkeypatch.setattr(versions.metadata, "distributions", boom)
    assert versions.python_package("numpy")["version"] is None


def test_a_resolver_that_fails_gives_unknown():
    def boom():
        raise OSError("/var/lib/mzdesk/secret is unreadable")
    entry = versions._guard("engine", boom)
    assert entry == {"version": None, "kind": "engine", "detail": "unknown"}, "the error text is not passed on"


def test_a_resolver_cannot_leak_a_path():
    entry = versions._guard("data", lambda: {"version": "1", "kind": "data", "detail": "built in /var/lib/mzdesk/data/libraries/GNPS on 2026-10-30"})
    assert "/var/lib" not in json.dumps(entry) and "GNPS on 2026-10-30" in entry["detail"]


def test_nothing_installed_gives_null_versions(nothing_installed):
    components = versions.environment()["components"]
    for name in ("MS-FINDER", "FRIGID", "CFM-ID", "MSRB", "MetFrag", "Mono", "Java", "NEIMS", "MolGpKa",
                 "MassBank", "GNPS", "MassBank/GNPS", "CCSbase", "RepoRT"):
        assert components[name]["version"] is None, name
        assert components[name]["detail"].startswith("not installed"), (name, components[name]["detail"])
    assert components["CFM-ID"]["detail"] == "not installed (Docker not found)"
    assert components["MetFrag"]["detail"] == "not installed (Java not found)"
    assert components["mzDesk"]["version"] == __version__ and components["numpy"]["version"], "what is installed is still found"


def test_a_library_index_gives_its_release_and_size_but_not_its_folder(nothing_installed, monkeypatch, tmp_path):
    from mzdesk.ms2 import library_index
    folder = tmp_path / "libraries"
    (folder / "MassBank").mkdir(parents=True)
    (folder / "MassBank" / "meta.json").write_text(json.dumps({"format": library_index.FORMAT, "library": "MassBank", "version": "2026.03", "date": "2026-04-15",
                                                                 "built": "2026-10-30 01:29:23", "n_spectra": 126507}))
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(folder))
    env = versions.environment()
    massbank = env["components"]["MassBank"]
    assert massbank == {"version": "2026.03", "kind": "data", "detail": "126,507 spectra, built 2026-10-30, source dated 2026-04-15"}
    assert env["components"]["GNPS"]["version"] is None
    both = env["components"]["MassBank/GNPS"]
    assert both["version"] == "MassBank 2026.03" and "GNPS: not installed" in both["detail"]
    assert str(tmp_path) not in json.dumps(env)


def test_an_old_library_index_is_reported_as_not_installed(nothing_installed, monkeypatch, tmp_path):
    (tmp_path / "libs" / "GNPS").mkdir(parents=True)
    (tmp_path / "libs" / "GNPS" / "meta.json").write_text(json.dumps({"format": 1, "version": "x", "n_spectra": 3}))
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "libs"))
    assert versions.environment()["components"]["GNPS"] == {"version": None, "kind": "data", "detail": "not installed (the index has an old format)"}


def test_model_files_give_the_model_version_and_file_name_only(nothing_installed, monkeypatch, tmp_path):
    neims, molgpka = tmp_path / "neims_weights.npz", tmp_path / "molgpka_weights.npz"
    neims.write_bytes(b"x")
    molgpka.write_bytes(b"x")
    monkeypatch.setenv("MZDESK_NEIMS_WEIGHTS", str(neims))
    monkeypatch.setenv("MZDESK_MOLGPKA_WEIGHTS", str(molgpka))
    env = versions.environment()
    from mzdesk.platform_fit import pka_setup
    assert env["components"]["NEIMS"]["version"] == "model.ckpt-100000" and "(neims_weights.npz)" in env["components"]["NEIMS"]["detail"]
    assert env["components"]["MolGpKa"]["version"] == pka_setup.COMMIT[:7] and "(molgpka_weights.npz)" in env["components"]["MolGpKa"]["detail"]
    assert str(tmp_path) not in json.dumps(env)


# ---------------------------------------------------------------- the three fixed commands
def fake_commands(monkeypatch, answers):
    """Replace versions._run: answers[program name] = (exit code, output); every call is recorded."""
    calls = []

    def run(cmd, timeout):
        calls.append((list(cmd), timeout))
        return answers.get(Path(cmd[0]).name)
    monkeypatch.setattr(versions, "_run", run)
    monkeypatch.setattr(versions.shutil, "which", lambda name, *a, **k: f"/fake/{name}" if name in ("docker", "java") else None)
    return calls


def test_docker_mono_and_java_are_asked_with_fixed_commands_and_a_time_limit(monkeypatch, tmp_path):
    mono = tmp_path / "mono"
    mono.write_text("#!/bin/sh\n")
    mono.chmod(0o755)
    monkeypatch.setenv("MZDESK_MONO", str(mono))
    from mzdesk.ms2 import msfinder_engine
    monkeypatch.setattr(msfinder_engine, "find_mono", lambda: str(mono))     # on Windows MS-FINDER runs without Mono: ask it anyway here
    calls = fake_commands(monkeypatch, {
        "docker": (0, "sha256:144256468e096765d6236f07f4c304e1d04d20dcf891aa2008cfec2fa9bbe1f8\n"),
        "java": (0, 'openjdk version "21.0.13.0.3" 2026-10-24\nOpenJDK Runtime Environment (build 21.0.13.0.3)\n'),
        "mono": (0, "Mono JIT compiler version 6.12.0.182 (tarball Wed Mar 6 2024)\nCopyright\n")})
    components = versions.environment()["components"]
    commands = sorted(c[0][0] if Path(c[0][0]).name != "mono" else "mono" for c in calls)
    assert commands == ["docker", "java", "mono"], "each command runs once per answer, however many components read it"
    for cmd, timeout in calls:
        assert timeout <= 5
        if cmd[0] == "docker":
            assert cmd == ["docker", "image", "inspect", "--format", "{{.Id}}", cfmid_engine.IMAGE]
        elif cmd[0] == "java":
            assert cmd == ["java", "-version"]
        else:
            assert cmd == [str(mono), "--version"]
    assert components["CFM-ID"] == {"version": "4.4.10", "kind": "engine", "detail": f"Docker image {cfmid_engine.IMAGE} (image id 144256468e09)"}
    assert components["MSRB"]["version"] == "1.1.15"
    assert components["Java"]["version"] == "21.0.13.0.3" and components["Mono"]["version"] == "6.12.0.182"
    assert str(mono) not in json.dumps(components), "the path of Mono stays on the server"


def test_docker_without_the_image_or_the_daemon_gives_null(monkeypatch):
    fake_commands(monkeypatch, {"docker": (1, "Error: No such image: wishartlab/cfmid:x")})
    assert versions.environment()["components"]["CFM-ID"] == {"version": None, "kind": "engine", "detail": "not installed (Docker image not pulled)"}
    versions.reset()
    versions._docker_image.cache_clear()
    fake_commands(monkeypatch, {"docker": (1, "Cannot connect to the Docker daemon at unix:///var/run/docker.sock")})
    entry = versions.environment()["components"]["MSRB"]
    assert entry == {"version": None, "kind": "engine", "detail": "not installed (Docker is not running)"}
    versions.reset()
    versions._docker_image.cache_clear()
    fake_commands(monkeypatch, {})
    assert versions.environment()["components"]["CFM-ID"]["detail"] == "not installed (Docker did not answer)"


def test_run_survives_a_missing_program_and_a_timeout():
    assert versions._run(["no-such-program-xyz"], 1) is None
    assert versions._run([sys.executable, "-c", "import time; time.sleep(5)"], 0.3) is None
    code, text = versions._run([sys.executable, "-c", "print('hi')"], 10)
    assert code == 0 and text.strip() == "hi"


# ---------------------------------------------------------------- the answer is kept
def test_the_answer_is_cached_and_looked_up_again_after_the_time_limit(monkeypatch):
    calls = fake_commands(monkeypatch, {"java": (0, 'java version "1.8.0_392"')})
    first = versions.environment()
    n = len(calls)
    assert n >= 1
    assert versions.environment() == first and len(calls) == n, "the second call does not look again"
    first["components"]["numpy"]["version"] = "tampered"
    assert versions.environment()["components"]["numpy"]["version"] != "tampered", "a caller gets its own copy"
    monkeypatch.setattr(versions, "TTL_S", -1)
    versions.environment()
    assert len(calls) > n, "after the time limit the next call looks again"


# ---------------------------------------------------------------- what the metadata says, and an upgrade that left old metadata behind
def make_dist(site, name, version, module_source):
    """<site>/<name>/__init__.py and a <name>-<version>.dist-info whose RECORD lists that file with the hash of `module_source` (as pip writes it)."""
    import base64
    import hashlib
    folder = site / f"{name}-{version}.dist-info"
    folder.mkdir(parents=True)
    (folder / "METADATA").write_text(f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n")
    digest = base64.urlsafe_b64encode(hashlib.sha256(module_source.encode()).digest()).rstrip(b"=").decode()
    (folder / "RECORD").write_text(f"{name}/__init__.py,sha256={digest},{len(module_source)}\n")


def test_a_package_that_forgot_to_update_its_own_version_is_read_from_the_metadata(monkeypatch):
    # msbuddy 0.3.12 still says __version__ = "0.3.9": the version is the installed distribution's, loaded or not
    monkeypatch.setitem(sys.modules, "fake_stale_xyz", types.SimpleNamespace(__version__="0.3.9", __file__=__file__))
    monkeypatch.setattr(versions.metadata, "distributions", lambda **kw: iter([types.SimpleNamespace(version="0.3.12", files=None)]))
    assert versions._dist_version("fake-stale-xyz", module="fake_stale_xyz") == "0.3.12"
    assert versions.python_package("fake-stale-xyz", module="fake_stale_xyz")["version"] == "0.3.12"


def test_the_real_buddy_version_is_the_distribution_whether_or_not_msbuddy_is_imported():
    listed = {d.version for d in versions.metadata.distributions(name="msbuddy")}
    before = versions.python_package("msbuddy", module="msbuddy")["version"]
    import msbuddy                                # msbuddy.__version__ may be out of date; the answer must not change
    assert before in listed and versions.python_package("msbuddy", module="msbuddy")["version"] == before


def test_two_metadata_folders_of_one_package_give_the_release_that_wrote_the_imported_file(monkeypatch, tmp_path):
    import importlib
    old_source, new_source = "VERSION = 1\n", "VERSION = 2  # newer\n"
    make_dist(tmp_path, "fakepkg_dup", "1.0", old_source)
    make_dist(tmp_path, "fakepkg_dup", "2.0", new_source)
    (tmp_path / "fakepkg_dup").mkdir()
    (tmp_path / "fakepkg_dup" / "__init__.py").write_bytes(new_source.encode())
    monkeypatch.syspath_prepend(str(tmp_path))
    importlib.invalidate_caches()
    for loaded in (False, True):                  # found with or without importing it
        if loaded:
            monkeypatch.delitem(sys.modules, "fakepkg_dup", raising=False)
            monkeypatch.setitem(sys.modules, "fakepkg_dup", importlib.import_module("fakepkg_dup"))
        assert versions._dist_version("fakepkg-dup", module="fakepkg_dup") == "2.0"
    # the file on disk is the OLD release (an upgrade that did not finish): the highest version listed must not win
    (tmp_path / "fakepkg_dup" / "__init__.py").write_bytes(old_source.encode())
    assert versions._dist_version("fakepkg-dup", module="fakepkg_dup") == "1.0"
    # neither record fits the file (edited by hand): the highest version listed
    (tmp_path / "fakepkg_dup" / "__init__.py").write_text("VERSION = 3\n")
    assert versions._dist_version("fakepkg-dup", module="fakepkg_dup") == "2.0"
    # no module name given: the highest version listed
    (tmp_path / "fakepkg_dup" / "__init__.py").write_bytes(old_source.encode())
    assert versions._dist_version("fakepkg-dup") == "2.0"


def test_scikit_learn_with_two_metadata_folders_is_the_one_imported():
    import sklearn
    from importlib.metadata import distributions
    listed = {d.version for d in distributions(name="scikit-learn")}
    got = versions._dist_version("scikit-learn", module="sklearn")
    assert got in listed
    if len(listed) > 1:                           # this machine: 1.7.1 and 1.8.0 side by side
        assert got == sklearn.__version__


# ---------------------------------------------------------------- nothing about the server leaves
def test_the_answer_holds_no_folder_user_or_secret():
    text = json.dumps(client.get("/api/versions").json())
    assert not ABSOLUTE_PATH.search(text), ABSOLUTE_PATH.search(text)
    assert str(Path.home()) not in text and os.getcwd() not in text
    for var, value in os.environ.items():
        if re.search(r"KEY|SECRET|PASSWORD|TOKEN|HASH|SALT", var) and len(value) >= 6:
            assert value not in text, var
