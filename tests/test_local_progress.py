"""Progress lines of the install steps (mzdesk/local/progress.py) and how a job row follows them."""
import sys

from mzdesk.local import jobs, progress


def test_count_line_gives_the_fraction():
    assert progress.parse("  412 / 2300 MB") == {"got_mb": 412.0, "total_mb": 2300.0, "progress": 0.1791}
    assert progress.parse("1 / 2 GB")["total_mb"] == 2048.0


def test_a_new_download_resets_the_numbers_and_names_the_file():
    got = progress.parse("Downloading https://x.org/a/b/MSFINDER.v3.73.zip?token=1")
    assert got == {"step": "Downloading MSFINDER.v3.73.zip", "progress": None, "got_mb": None, "total_mb": None}


def test_phases_and_unknown_lines():
    assert progress.parse("Extracting DLM.ckpt from DLM.tar.gz")["step"] == "Unpacking"
    assert progress.parse("Compiling /x/MsfinderConsoleApp.exe")["step"] == "Compiling"
    assert progress.parse("Installing the pinned wheels (about 2 GB)")["step"] == "Installing Python packages"
    assert progress.parse("nothing to see") is None
    assert progress.parse("  5 / 0 MB") is None


def test_a_running_row_follows_the_lines(monkeypatch):
    job = jobs.Job([], None)
    job.rows = [{"id": "a", "title": "A", "status": "running", "detail": "", "step": "", "progress": None, "got_mb": None, "total_mb": None},
                {"id": "b", "title": "B", "status": "waiting", "detail": "", "step": "", "progress": None, "got_mb": None, "total_mb": None}]
    job._write("Downloading https://h/f.zip")
    job._write("  50 / 200 MB")
    row = job.snapshot()["components"][0]
    assert (row["step"], row["progress"], row["got_mb"]) == ("Downloading f.zip", 0.25, 50.0)
    assert job.snapshot()["components"][1]["progress"] is None          # the waiting row is untouched
    job._set("a", "done", "installed")
    assert job.snapshot()["components"][0]["progress"] is None and job.snapshot()["components"][0]["step"] == ""
