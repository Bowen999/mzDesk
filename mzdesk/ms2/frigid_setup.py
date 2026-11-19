"""FRIGID-base 安装与检查，一次性运行：python -m mzdesk.ms2.frigid_setup [--check] [--verify] [--force] [--dir 文件夹]

De novo structure 页用 FRIGID-base（Coley 组，Apache-2.0）在本机生成候选结构。这个命令把它装到 ~/.mzdesk/tools/frigid/（MZDESK_FRIGID_DIR 可改；~/tools/frigid 已存在时用它），共约 3 GB 下载：
  1. env/        conda 前缀环境 Python 3.10（conda-forge，连同预编译的 antlr 运行库）+ 钉死版本的轮子（mzdesk/lib/frigid/requirements.txt，含 PyTorch 2.6）；
  2. FRIGID/     github.com/coleygroup/FRIGID 的克隆，固定在提交 9a3e623，打两个补丁：mps_patch.diff（设备选择 cuda → mps → cpu）和 SAFE 包的导入补丁
                 （FRIGID 自己的 env/fix_safe_imports.sh 的改动，这里用 Python 重做，因为原脚本用 GNU sed）；
  3. checkpoints/mist_msg_aug.ckpt   Zenodo 23303289（MIT，289 MB，校验 md5）；
     checkpoints/DLM.ckpt            Zenodo 23543494 的 tar.gz（Apache-2.0，2.6 GB，校验 md5），只取出 DLM.ckpt，随后删掉 tar；
  4. hf_cache/   SAFE 分词器 datamol-io/safe-gpt 的 tokenizer.json（46 KB，固定到一个提交），运行时离线读取。
驱动脚本（mzdesk/ms2/frigid_driver/）和 SAFE 长度模型（mzdesk/lib/frigid/，我们在 MIT 许可的 MassSpecGym 标签上训练）随 mzDesk 一起发布，不在安装文件夹里。
每一步都可重复运行：已经就位（大小 / md5 / 提交都对）就跳过。只有安装时联网；--check 只报告状态，--verify 再算两个检查点的 md5。
引用：FRIGID，Coley 组 https://github.com/coleygroup/FRIGID ；MIST，Goldman et al.，github.com/samgoldman97/mist（MIT）。
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

from mzdesk.ms2 import frigid_engine as fe

REPO_URL = "https://github.com/coleygroup/FRIGID"
PYTHON_VERSION = "3.10"
ANTLR_CONDA = "antlr-python-runtime=4.9.3"          # the PyPI package is source-only; conda-forge has a prebuilt one

DLM_URL = "https://zenodo.org/api/records/23543494/files/frigid_pretrained_checkpoints.tar.gz/content"
DLM_TAR = "frigid_pretrained_checkpoints.tar.gz"
DLM_TAR_SIZE = 2_614_091_122
DLM_TAR_MD5 = "5a663a3000664ded0cc62fc5a4ca39e0"
DLM_MD5 = "717f93448661d9b63c4ceb4f7bcdfab0"
MIST_URL = "https://zenodo.org/api/records/23303289/files/mist_msg_aug.ckpt/content"
MIST_MD5 = "3e92179240e3043733627d720f4b7b86"
TOKENIZER_SIZE = 46_319
SAFE_INIT_LAST_LINE = "from .converter import SAFEConverter, decode, encode"


class SetupError(RuntimeError):
    pass


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify_file(path, size, md5):
    """Raise SetupError unless the file has the pinned size and md5."""
    path = Path(path)
    if not path.is_file():
        raise SetupError(f"{path} is missing.")
    if path.stat().st_size != size:
        raise SetupError(f"{path} has {path.stat().st_size} bytes, expected {size}.")
    digest = md5_of(path)
    if digest != md5:
        raise SetupError(f"{path} has md5 {digest}, expected {md5}.")


def file_ok(path, size, md5):
    try:
        verify_file(path, size, md5)
    except SetupError:
        return False
    return True


def _run(cmd, **kw):
    """subprocess.run with text output; a missing program becomes a SetupError."""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, **kw)
    except OSError as exc:
        raise SetupError(f"Could not run {cmd[0]}: {exc}") from None


def _env():
    return {**os.environ, "PYTHONNOUSERSITE": "1"}


# ---------------------------------------------------------------- download
def download(url, dest, expect_size=None):
    """Download ``url`` to ``dest`` through ``dest.part`` (an interrupted download resumes when the server supports Range requests)."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    offset = part.stat().st_size if part.exists() else 0
    print(f"Downloading {url}\n  -> {dest}" + (f" ({expect_size / 1e6:.0f} MB, this can take a while)" if expect_size else ""))
    headers = {"User-Agent": "mzdesk-frigid-setup"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=60) as resp:
            resume = bool(offset) and getattr(resp, "status", None) == 206
            done, report = (offset if resume else 0), 10 << 20
            total = expect_size or (done + int(resp.headers.get("Content-Length") or 0))
            with open(part, "ab" if resume else "wb") as out:
                while block := resp.read(1 << 20):
                    out.write(block)
                    done += len(block)
                    if done >= report:
                        print(f"  {done >> 20} / {total >> 20} MB" if total > done else f"  {done >> 20} MB", flush=True)
                        report = done + (10 << 20)
    except OSError as exc:
        raise SetupError(f"Download failed ({exc}). Run the command again: a partial download is resumed.") from None
    part.replace(dest)


def fetch_file(url, dest, size, md5, force=False):
    """Download and verify one file; skip when it is already there and correct. Returns True when a download happened."""
    dest = Path(dest)
    if dest.exists() and not force:
        if file_ok(dest, size, md5):
            return False
        raise SetupError(f"{dest} exists but its size or md5 is wrong. Delete it (or use --force) and run again.")
    if force and dest.exists():
        dest.unlink()
    download(url, dest, size)
    try:
        verify_file(dest, size, md5)
    except SetupError:
        dest.replace(dest.with_name(dest.name + ".bad"))
        raise
    return True


def extract_member(tar_path, name, dest, size, md5):
    """Copy the one member whose file name is ``name`` out of a (gzipped) tar into ``dest``, verifying size and md5. Nothing else is extracted."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    h, n = hashlib.md5(), 0
    with tarfile.open(tar_path, "r|*") as tf:                      # streaming: the 2.6 GB archive is read once, never held in memory
        for member in tf:
            if member.isfile() and Path(member.name).name == name:
                with tf.extractfile(member) as src, open(part, "wb") as out:
                    for block in iter(lambda: src.read(1 << 20), b""):
                        out.write(block)
                        h.update(block)
                        n += len(block)
                break
        else:
            raise SetupError(f"{name} is not in {tar_path}.")
    if n != size or h.hexdigest() != md5:
        part.unlink()
        raise SetupError(f"{name} from {tar_path} has {n} bytes and md5 {h.hexdigest()}, expected {size} bytes and {md5}.")
    part.replace(dest)


# ---------------------------------------------------------------- steps
def conda_candidates():
    """Where a conda / mamba may be, besides PATH and CONDA_EXE (the folders differ on Windows: Scripts\\conda.exe)."""
    home = Path.home()
    if sys.platform == "win32":
        return [str(home / base / "Scripts" / "conda.exe") for base in ("miniforge3", "miniconda3", "anaconda3", "mambaforge", "micromamba")] + [
            str(home / "AppData" / "Local" / base / "Scripts" / "conda.exe") for base in ("miniforge3", "miniconda3", "anaconda3")] + [
            str(Path(os.environ.get("ProgramData", r"C:\ProgramData")) / base / "Scripts" / "conda.exe") for base in ("miniforge3", "miniconda3", "anaconda3")]
    return ["/opt/anaconda3/bin/conda", str(home / "miniforge3/bin/conda"), str(home / "miniconda3/bin/conda"),
            str(home / "anaconda3/bin/conda"), "/opt/homebrew/Caskroom/miniforge/base/bin/conda"]


def find_conda():
    for cand in (os.environ.get("CONDA_EXE"), shutil.which("conda"), shutil.which("mamba"), shutil.which("micromamba"), *conda_candidates()):
        if cand and Path(cand).is_file() and os.access(cand, os.X_OK):
            return cand
    return None


def env_create_command(conda, env_dir):
    return [conda, "create", "-y", "-p", str(env_dir), "--override-channels", "-c", "conda-forge", f"python={PYTHON_VERSION}", ANTLR_CONDA]


def pip_install_command(py, requirements=None):
    return [str(py), "-m", "pip", "install", "--only-binary=:all:", "-r", str(requirements or fe.REQUIREMENTS)]


def ensure_env(install):
    """Create <install>/env (Python 3.10 + antlr from conda-forge) when it is missing; install the pinned wheels either way (pip skips what is there)."""
    py = fe.env_python(install)
    if not py.is_file():
        conda = find_conda()
        if not conda:
            raise SetupError("conda was not found. Install Miniforge or Anaconda (or put conda on PATH / set CONDA_EXE), then run this command again. "
                             f"The environment needs Python {PYTHON_VERSION}.")
        print(f"Creating the Python environment {install / 'env'}")
        r = _run(env_create_command(conda, install / "env"))
        if r.returncode != 0 or not py.is_file():
            raise SetupError("conda create failed:\n" + (r.stdout + r.stderr)[-3000:])
    print("Installing the pinned wheels (about 2 GB incl. PyTorch; already installed ones are skipped)")
    r = subprocess.run(pip_install_command(py), env=_env())
    if r.returncode != 0:
        raise SetupError("pip install failed (see the output above). Run the command again.")
    return py


def clone_frigid(install):
    """Clone FRIGID at the pinned commit; an existing clone must already be at that commit (it is never moved)."""
    clone = install / "FRIGID"
    if (clone / ".git").exists():
        head = _run(["git", "-C", str(clone), "rev-parse", "HEAD"]).stdout.strip()
        if head != fe.FRIGID_COMMIT:
            raise SetupError(f"{clone} is at {head or 'an unknown commit'}, mzDesk is verified with {fe.FRIGID_COMMIT}. "
                             "Move or delete the folder and run again.")
        return clone, False
    if clone.exists() and any(clone.iterdir()):
        raise SetupError(f"{clone} exists and is not a git clone.")
    for cmd in (["git", "clone", "--quiet", "--config", "core.autocrlf=false", REPO_URL, str(clone)], ["git", "-C", str(clone), "checkout", "--quiet", "--detach", fe.FRIGID_COMMIT]):
        r = _run(cmd)
        if r.returncode != 0:
            raise SetupError(f"{' '.join(cmd[:3])} failed:\n" + (r.stdout + r.stderr)[-2000:])
    return clone, True


def apply_mps_patch(clone, patch=None):
    """Apply mps_patch.diff to the clone unless it is already applied. Returns "applied" | "already applied"."""
    patch = str(patch or fe.MPS_PATCH)
    if _run(["git", "-C", str(clone), "apply", "--reverse", "--check", patch]).returncode == 0:
        return "already applied"
    r = _run(["git", "-C", str(clone), "apply", "--check", patch])
    if r.returncode != 0:
        raise SetupError("The device patch does not apply to the clone (is it at the pinned commit with no other edits?):\n" + (r.stdout + r.stderr)[-1500:])
    r = _run(["git", "-C", str(clone), "apply", patch])
    if r.returncode != 0:
        raise SetupError("git apply failed:\n" + (r.stdout + r.stderr)[-1500:])
    return "applied"


def apply_safe_patch(install):
    """FRIGID's env/fix_safe_imports.sh edit to <env>/.../safe/__init__.py: comment out every line, append one import. Idempotent; keeps __init__.py.orig."""
    init = fe.safe_init(install)
    if init is None:
        raise SetupError(f"The safe package is not installed in {install / 'env'}.")
    if fe.safe_patch_applied(init):
        return "already applied"
    orig = init.with_name("__init__.py.orig")
    if not orig.exists():
        shutil.copy2(init, orig)
    with open(init, "w", encoding="utf-8", newline="\n") as f:
        f.write("".join("# " + ln + "\n" for ln in init.with_name("__init__.py.orig").read_text(encoding="utf-8").splitlines()) + SAFE_INIT_LAST_LINE + "\n")
    return "applied"


def fetch_checkpoints(install, force=False):
    ckpt = install / "checkpoints"
    if fetch_file(MIST_URL, ckpt / "mist_msg_aug.ckpt", fe.MIST_SIZE, MIST_MD5, force):
        print("mist_msg_aug.ckpt downloaded and verified")
    dlm = ckpt / "DLM.ckpt"
    if not force and file_ok(dlm, fe.DLM_SIZE, DLM_MD5):
        return
    tar = install / "downloads" / DLM_TAR
    fetch_file(DLM_URL, tar, DLM_TAR_SIZE, DLM_TAR_MD5, force)
    print(f"Extracting DLM.ckpt from {tar.name} (the other checkpoints in it are not needed)")
    extract_member(tar, "DLM.ckpt", dlm, fe.DLM_SIZE, DLM_MD5)
    tar.unlink()
    print("DLM.ckpt verified; the tar.gz is deleted")


def write_main_ref(hf_cache, revision=None):
    """hf_cache/hub/models--datamol-io--safe-gpt/refs/main = the pinned revision (a download by commit hash does not write it; the offline lookup of
    revision "main" needs it)."""
    ref = Path(hf_cache) / "hub" / "models--datamol-io--safe-gpt" / "refs" / "main"
    if not ref.exists():
        ref.parent.mkdir(parents=True, exist_ok=True)
        ref.write_text(revision or fe.TOKENIZER_REVISION, encoding="utf-8")
        return True
    return False


def fetch_tokenizer(install):
    """Cache the 46 KB SAFE tokenizer (datamol-io/safe-gpt, pinned revision) in hf_cache with the env's huggingface_hub. Skipped when it is there."""
    cached = fe.tokenizer_file(install)
    if cached and cached.stat().st_size == TOKENIZER_SIZE:
        write_main_ref(install / "hf_cache")
        return False
    code = f"from huggingface_hub import hf_hub_download as d; print(d({fe.TOKENIZER_REPO!r}, 'tokenizer.json', revision={fe.TOKENIZER_REVISION!r}))"
    env = {**_env(), "HF_HOME": str(install / "hf_cache")}
    env.pop("HF_HUB_OFFLINE", None)
    r = _run([str(fe.env_python(install)), "-I", "-c", code], env=env)
    cached = fe.tokenizer_file(install)
    if r.returncode != 0 or not cached or cached.stat().st_size != TOKENIZER_SIZE:
        raise SetupError("Could not fetch the SAFE tokenizer from the Hugging Face Hub:\n" + (r.stdout + r.stderr)[-1500:])
    write_main_ref(install / "hf_cache")
    return True


# ---------------------------------------------------------------- status and setup
def check(verify=False):
    """Print what is in place; returns 0 when De novo structure can run, else 1. Starts nothing; --verify also computes the checkpoints' md5."""
    s = fe.status()
    print("FRIGID status (mzDesk De novo structure)")
    print(f"  Install folder     : {s['dir']}  (MZDESK_FRIGID_DIR changes it)")
    for key, label, missing in (("python", "Python environment", "env/bin/python"), ("clone", "FRIGID code", f"FRIGID/ at {fe.FRIGID_COMMIT[:7]}"),
                                ("dlm", "DLM checkpoint", "checkpoints/DLM.ckpt"), ("mist", "MIST checkpoint", "checkpoints/mist_msg_aug.ckpt"),
                                ("tokenizer", "SAFE tokenizer", "hf_cache/")):
        print(f"  {label:<19}: {s[key] if s[key] else 'NOT FOUND (' + missing + ')'}")
    print(f"  Device patch       : {'applied' if s['mps_patch'] else 'NOT APPLIED (FRIGID/src/dlm/sampler.py)'}")
    print(f"  SAFE import patch  : {'applied' if s['safe_patch'] else 'NOT APPLIED (env: safe/__init__.py)'}")
    print(f"  mzDesk driver      : {s['driver'] if s['driver'] else 'NOT FOUND (reinstall mzDesk)'}")
    print(f"  SAFE length model  : {s['length_model'] if s['length_model'] else 'NOT FOUND (reinstall mzDesk)'}")
    ready = all(s[k] for k in ("python", "clone", "mps_patch", "safe_patch", "dlm", "mist", "tokenizer", "driver", "length_model"))
    if verify:
        for key, md5 in (("dlm", DLM_MD5), ("mist", MIST_MD5)):
            if s[key]:
                ok = md5_of(s[key]) == md5
                print(f"  md5 {s[key].name:<19}: {'ok' if ok else 'WRONG (delete the file and run the setup again)'}")
                ready = ready and ok
    print("  Ready              : " + ("yes" if ready else f"no. Run: {fe.SETUP_COMMAND}"))
    return 0 if ready else 1


def setup(install, force=False):
    print("Installing FRIGID-base for De novo structure (FRIGID code Apache-2.0, DLM checkpoint Apache-2.0, MIST checkpoint MIT).")
    install.mkdir(parents=True, exist_ok=True)
    ensure_env(install)
    clone, cloned = clone_frigid(install)
    print(f"FRIGID clone: {clone} ({'cloned' if cloned else 'present'}, commit {fe.FRIGID_COMMIT[:7]})")
    print(f"Device patch: {apply_mps_patch(clone)}")
    print(f"SAFE import patch: {apply_safe_patch(install)}")
    fetch_checkpoints(install, force)
    print("SAFE tokenizer: " + ("downloaded" if fetch_tokenizer(install) else "present"))
    return check()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Install or check FRIGID-base for the De novo structure page.")
    ap.add_argument("--check", action="store_true", help="only report what is installed (no download, nothing is started)")
    ap.add_argument("--verify", action="store_true", help="with --check: also compute the md5 of the two checkpoints (a few seconds)")
    ap.add_argument("--force", action="store_true", help="download the checkpoints again even if present")
    ap.add_argument("--dir", help="install folder (default: MZDESK_FRIGID_DIR, ~/tools/frigid if it exists, else ~/.mzdesk/tools/frigid)")
    args = ap.parse_args(argv)
    if args.dir:
        os.environ["MZDESK_FRIGID_DIR"] = args.dir
    if args.check or args.verify:
        return check(verify=args.verify)
    try:
        return setup(fe.install_dir(), args.force)
    except SetupError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
