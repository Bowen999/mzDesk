"""What a log line of an install step says about progress. Standard library only.

The download helpers of the setup commands print one line per file when it starts ("Downloading <url>") and, while it runs, `  412 / 2300 MB`.
parse() turns a line into a change of the row's record: a new step name, or how far the file is. Lines it does not know change nothing.

    parse("  412 / 2300 MB")           -> {"got_mb": 412.0, "total_mb": 2300.0, "progress": 0.179}
    parse("Downloading https://x/y.zip") -> {"step": "Downloading y.zip", "progress": None, "got_mb": None, "total_mb": None}
"""
import re

_UNIT_MB = {"KB": 1 / 1024, "MB": 1.0, "GB": 1024.0}
_COUNT = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)\s*(KB|MB|GB)\s*$")
_DOWNLOAD = re.compile(r"^Downloading\s+(\S+)")
_STEPS = (                                            # first match wins; the text is what the dialog shows next to the bar
    (re.compile(r"extract|unzip|unpack", re.I), "Unpacking"),
    (re.compile(r"compil", re.I), "Compiling"),
    (re.compile(r"building the .*index|^index\b|rebuilding", re.I), "Building the search index"),
    (re.compile(r"python environment|pinned wheels|pip install|installing .*(wheels|packages)", re.I), "Installing Python packages"),
    (re.compile(r"cloning|checkout|sparse", re.I), "Fetching data"),
    (re.compile(r"tensor|convert|weights", re.I), "Converting model weights"),
    (re.compile(r"verified|checking|sha-?256|md5", re.I), "Verifying"),
    (re.compile(r"docker (pull|image)|pulling", re.I), "Pulling the Docker image"),
    (re.compile(r"training|trains? the", re.I), "Training the model"),
)


def _file_name(url):
    name = url.split("?", 1)[0].rstrip("/").rsplit("/", 1)[-1]
    return name[:60] if name else url[:60]


def parse(line):
    """-> dict of the row fields this line changes, or None."""
    m = _COUNT.match(line)
    if m:
        unit = _UNIT_MB[m.group(3)]
        got, total = float(m.group(1)) * unit, float(m.group(2)) * unit
        if total > 0:
            return {"got_mb": round(got, 1), "total_mb": round(total, 1), "progress": round(min(1.0, got / total), 4)}
        return None
    d = _DOWNLOAD.match(line)
    if d:
        return {"step": f"Downloading {_file_name(d.group(1))}", "progress": None, "got_mb": None, "total_mb": None}
    text = line.strip()
    for pattern, step in _STEPS:
        if pattern.search(text):
            return {"step": step, "progress": None, "got_mb": None, "total_mb": None}
    return None
