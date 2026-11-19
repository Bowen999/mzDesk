"""Download the files of one pinned GitHub commit without git: stream the commit's tar.gz from codeload.github.com and keep only the wanted files.

Used by the CCS and RepoRT setups (the data they need is a few files of a big repository). The commit is part of the URL, so the content is
exactly that commit's; the archive's top folder name (<repo>-<commit>) is checked as well. Standard library only.
"""
import re
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

MARKER = ".mzdesk-commit"          # written when all files of a commit are in place
_UA = {"User-Agent": "mzdesk-setup"}


class ArchiveError(RuntimeError):
    pass


def archive_url(repo_url, commit):
    """https://github.com/<owner>/<repo>[.git] + commit -> https://codeload.github.com/<owner>/<repo>/tar.gz/<commit>"""
    m = re.fullmatch(r"https://github\.com/([\w.-]+)/([\w.-]+?)(?:\.git)?/?", repo_url)
    if not m or not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ArchiveError(f"Not a GitHub repository URL or not a full commit hash: {repo_url} {commit}")
    return f"https://codeload.github.com/{m.group(1)}/{m.group(2)}/tar.gz/{commit}"


def pattern_regex(pattern):
    """One sparse-checkout style pattern ('/LICENSE', '/a/*/b_*.tsv'; '*' stays inside one folder level) -> compiled regex for a repo-relative path."""
    parts = re.split(r"(\*)", pattern.lstrip("/"))
    return re.compile("".join("[^/]*" if p == "*" else re.escape(p) for p in parts))


def has_commit(folder, commit):
    """True when ``folder`` holds the files of ``commit``: our marker, or a git checkout that is on it (an install made with the old sparse clone)."""
    folder = Path(folder)
    try:
        if (folder / MARKER).read_text(encoding="utf-8").strip() == commit:
            return True
    except OSError:
        pass
    try:
        return (folder / ".git" / "HEAD").read_text(encoding="utf-8").strip() == commit
    except OSError:
        return False


def fetch(repo_url, commit, dest, patterns, log=print, timeout=60):
    """Stream the archive of ``commit`` and write the files that match ``patterns`` under ``dest`` (same relative paths as in the repository).

    Returns the number of files written; raises ArchiveError on a network problem, a damaged archive or when nothing matched."""
    url = archive_url(repo_url, commit)
    regexes = [pattern_regex(p) for p in patterns]
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    root = dest.resolve()
    written, seen, report = 0, 0, 25 << 20
    log(f"Downloading {url}")
    try:
        resp = urllib.request.urlopen(urllib.request.Request(url, headers=_UA), timeout=timeout)
    except (urllib.error.URLError, OSError) as exc:
        raise ArchiveError(f"Download failed ({exc}). Check the network connection and run again.") from None
    top = None
    try:
        with resp, tarfile.open(fileobj=resp, mode="r|gz") as tf:
            for member in tf:
                name = member.name
                if top is None:
                    top = name.split("/", 1)[0]
                    if not top.endswith(commit):
                        raise ArchiveError(f"The archive is for {top}, expected a folder that ends with {commit}.")
                rel = name.split("/", 1)[1] if "/" in name else ""
                if not member.isfile() or not rel or not any(r.fullmatch(rel) for r in regexes):
                    continue
                target = (root / rel).resolve()
                if root not in target.parents:
                    raise ArchiveError(f"Unsafe path in the archive: {name}")
                target.parent.mkdir(parents=True, exist_ok=True)
                part = target.with_name(target.name + ".part")
                with tf.extractfile(member) as src, open(part, "wb") as out:
                    for block in iter(lambda: src.read(1 << 20), b""):
                        out.write(block)
                part.replace(target)
                written += 1
                seen += member.size
                if seen >= report:
                    log(f"  {seen >> 20} MB kept")
                    report += 25 << 20
    except (tarfile.TarError, EOFError, OSError, urllib.error.URLError) as exc:
        raise ArchiveError(f"The download was interrupted or damaged ({exc}). Check the network connection and run again.") from None
    if not written:
        raise ArchiveError(f"No file of the archive matched {', '.join(patterns)}.")
    (dest / MARKER).write_text(commit + "\n", encoding="utf-8")
    return written
