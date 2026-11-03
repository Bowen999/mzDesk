"""MS-FINDER 安装与检查，一次性运行：python -m mzdesk.ms2.msfinder_setup [--check]

MS-FINDER v3.73（RIKEN，Tsugawa 实验室）只发布 Windows GUI 版。这里在本机用 Mono 跑它：
  1. 检查 Mono（没有就停下，给安装命令）；
  2. 下载官方 MSFINDER.v3.73.zip（约 325 MB，校验大小和 SHA-256）并解压到 ~/.mzdesk/tools/msfinder/MSFINDER_v3.73/（MZDESK_MSFINDER_DIR 可改；~/tools/msfinder/MSFINDER_v3.73 已存在时用它）；
  3. 取 tag MSFINDER-v3.73 里 MsfinderConsoleAppCore 的 9 个源文件（逐个校验固定的 SHA-256，一个字节都不改）；
  4. 用 Mono 的 csc 对着 zip 里的 DLL 编译出 MsfinderConsoleApp.exe，放进 MSFINDER_v3.73/（程序从 exe 所在目录找 Resources/）；
  5. 线程修复：StructureFinder.dll 里被并行打分线程共用的 static 字段 maxBondDissociationEnergy 标为 [ThreadStatic]
     （原因见 msfinder_engine 顶部；补丁源码 mzdesk/lib/msfinder/ThreadStaticPatch.cs，用 Mono 自带 GAC 里的 Mono.Cecil 编译运行）。
     只对 SHA-256 与发布版一致的文件打补丁；原文件留作 StructureFinder.dll.original，旁边写一个记录补丁结果 SHA-256 的标记文件。
只有步骤 2、3 联网；已装好的旧安装再运行一次只做第 5 步。--check 只报告状态，不下载、不启动 MS-FINDER。
许可：程序（assembly）CC BY 4.0，源码 LGPL-3.0，来自 github.com/systemsomicslab/MsdialWorkbench。引用：
Tsugawa H et al., Hydrogen rearrangement rules: computational MS/MS fragmentation and structure elucidation using MS-FINDER software,
Anal. Chem. 2016, 88, 7946-7958.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

from mzdesk.ms2 import msfinder_engine as mf

TAG = "MSFINDER-v3.73"
ZIP_URL = f"https://github.com/systemsomicslab/MsdialWorkbench/releases/download/{TAG}/MSFINDER.v3.73.zip"
ZIP_SIZE = 324_800_202
ZIP_SHA256 = "bf993637e282f0db2fb5a44b741ca8cffa361ff9bd88d1273715cc1b9cb25138"
ZIP_TOP = "MSFINDER v3.73/"                    # the folder inside the zip (with a space)
SOURCE_URL = f"https://raw.githubusercontent.com/systemsomicslab/MsdialWorkbench/{TAG}/src/MSFINDER/MsfinderConsoleAppCore/"
# the 9 files of src/MSFINDER/MsfinderConsoleAppCore at the tag, with their SHA-256 (copied from the files as fetched; nothing is patched)
SOURCES = {
    "MsfinderConsoleAppCore.csproj": "8fcc7a9e90d70db1e76b8aa74b7e83bb1a7ad654a2881b688803af5f9d2a4632",
    "Process/AnnotateProcess.cs": "1565267498428cbf9ead9fc81bd4675b0e08c726031e2726a25605460d26bad0",
    "Process/FragmentGenerateProcess.cs": "1d3df1ff4feef51342467f06870ec65b052671a1751e2f841641ac93177f1ba1",
    "Process/MainProcess.cs": "5e4c7d00247d1ce9006f7c1f78141687fa77d05e88eaaa75d8cf956bd956f43e",
    "Process/MsSearchProcess.cs": "d810852b59da47c2fff4c8d6445cfae53d73ce60a6faeda06b710492aa8b4828",
    "Process/PredictProcess.cs": "115f558c6785e3273b1cf33ae01249daa877061319ef38b75ee6012a99d74fdc",
    "Process/SpectrumCuration.cs": "ff0a9501423adf251612911cda5cf7724c978fce8b4e9a6d2a8af20c51cffb10",
    "Program.cs": "b390e407dc7000f90965895e2bf9ed82a2553f67b8f7bd280ac0fb54d3b47748",
    "Test/MergeMsdialMsfinderResult.cs": "5616651d16c7a450c35ad783e1eeb9b5f72f4e8972860db1712810c238b7436d",
}
DLLS = ["Common", "MolecularFormulaFinder", "MsfinderCommon", "StructureFinder", "ClassyfireApi", "PugRestApi", "Database", "DataProcessing",
        "Mathematics", "NCDK", "MessagePack"]


class SetupError(RuntimeError):
    pass


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def work_dir(install):
    """Where the zip and the console sources are kept: next to the install folder."""
    return install.parent


# ---------------------------------------------------------------- download
def _download(url, dest, expect_size=None):
    part = dest.with_name(dest.name + ".part")
    print(f"Downloading {url}\n  -> {dest}" + (f" ({expect_size / 1e6:.0f} MB, this can take a while)" if expect_size else ""))
    request = urllib.request.Request(url, headers={"User-Agent": "mzdesk-msfinder-setup"})
    with urllib.request.urlopen(request, timeout=60) as resp, open(part, "wb") as out:
        done, report = 0, 50 << 20
        while block := resp.read(1 << 20):
            out.write(block)
            done += len(block)
            if done >= report:
                print(f"  {done >> 20} MB", flush=True)
                report += 50 << 20
    part.rename(dest)


def get_zip(install, force=False):
    """Path of the verified release zip; downloads it when missing."""
    dest = work_dir(install) / "MSFINDER.v3.73.zip"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if force or not dest.exists():
        _download(ZIP_URL, dest, ZIP_SIZE)
    if dest.stat().st_size != ZIP_SIZE:
        raise SetupError(f"{dest} has {dest.stat().st_size} bytes, expected {ZIP_SIZE}. Delete it and run again.")
    digest = sha256_of(dest)
    if digest != ZIP_SHA256:
        raise SetupError(f"{dest} has sha256 {digest}, expected {ZIP_SHA256}. Delete it and run again.")
    return dest


def unzip(zip_path, install):
    """Extract the release into ``install`` (the folder inside the zip is called 'MSFINDER v3.73')."""
    install.mkdir(parents=True, exist_ok=True)
    root = install.resolve()
    with zipfile.ZipFile(zip_path) as z:
        for info in z.infolist():
            if not info.filename.startswith(ZIP_TOP) or info.filename == ZIP_TOP:
                continue
            target = (root / info.filename[len(ZIP_TOP):]).resolve()
            if root not in target.parents and target != root:
                raise SetupError(f"Unsafe path in the zip: {info.filename}")
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with z.open(info) as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst)


def verify_sources(src_dir, table=None):
    """Raise SetupError unless every source file exists with its pinned SHA-256."""
    for rel, expect in (table or SOURCES).items():
        path = Path(src_dir) / rel
        if not path.is_file():
            raise SetupError(f"Missing source file {path}.")
        digest = sha256_of(path)
        if digest != expect:
            raise SetupError(f"{path} has sha256 {digest}, expected {expect}. The source on GitHub changed or the file is damaged; "
                             "the exe is only built from the pinned files.")


def fetch_sources(src_dir, force=False):
    """Fetch the 9 console source files at the tag; each is verified before it is kept."""
    src_dir = Path(src_dir)
    for rel, expect in SOURCES.items():
        path = src_dir / rel
        if path.is_file() and not force and sha256_of(path) == expect:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part")
        _download(SOURCE_URL + rel, part)
        digest = sha256_of(part)
        if digest != expect:
            part.unlink()
            raise SetupError(f"{rel} fetched from {SOURCE_URL} has sha256 {digest}, expected {expect}.")
        part.replace(path)
    verify_sources(src_dir)


# ---------------------------------------------------------------- compile
def find_csc(mono):
    for cand in (shutil.which("csc"), str(Path(mono).resolve().parent / "csc")):
        if cand and Path(cand).is_file():
            return cand
    return None


def find_facades(mono):
    """Mono's lib/mono/4.5/Facades (netstandard.dll, System.Runtime.dll), found from the mono executable."""
    prefix = Path(mono).resolve().parent.parent
    for cand in (prefix / "lib" / "mono" / "4.5" / "Facades", Path("/usr/lib/mono/4.5/Facades")):
        if cand.is_dir():
            return cand
    return None


def compile_command(csc, facades, bin_dir, src_dir, out):
    """The csc command line (same as the one that built the verified exe): the DLLs of the release plus Mono's facades."""
    refs = [f"-r:{Path(bin_dir) / (d + '.dll')}" for d in DLLS]
    refs += [f"-r:{Path(facades) / 'netstandard.dll'}", f"-r:{Path(facades) / 'System.Runtime.dll'}", "-r:System.Xml.Linq.dll", "-r:System.Numerics.dll"]
    src = Path(src_dir)
    files = [src / "Program.cs", *sorted((src / "Process").glob("*.cs")), *sorted((src / "Test").glob("*.cs"))]
    return [csc, "-nologo", "-target:exe", "-platform:anycpu", f"-out:{out}", *refs, *map(str, files)]


def build_console(install, src_dir, out=None):
    """Compile MsfinderConsoleApp.exe from the pinned sources (into ``out``; default: the install folder). Returns the exe path."""
    mono = mf.find_mono()
    if not mono:
        raise SetupError("Mono not found.")
    csc, facades = find_csc(mono), find_facades(mono)
    if not csc or not facades:
        raise SetupError("Mono's C# compiler (csc) or its facade assemblies were not found; install the full Mono package "
                         "(for example `brew install mono`).")
    verify_sources(src_dir)
    for d in DLLS:
        if not (install / f"{d}.dll").is_file():
            raise SetupError(f"{install / (d + '.dll')} is missing; unzip the release first.")
    out = Path(out) if out else install / mf.EXE_NAME
    cmd = compile_command(csc, facades, install, src_dir, out)
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if proc.returncode != 0 or not out.is_file():
        raise SetupError("csc failed:\n" + ((proc.stdout or "") + (proc.stderr or ""))[-3000:])
    return out


# ---------------------------------------------------------------- thread fix of StructureFinder.dll
PATCH_SOURCE = mf.INI_PATH.parent / "ThreadStaticPatch.cs"
CECIL_VERSION = "0.11.1.0"              # in the GAC of Homebrew Mono 6.14 and of conda-forge Mono 6.12 (osx-64 and linux-64)


def find_cecil(mono):
    """Mono.Cecil.dll from the GAC of this Mono (CECIL_VERSION first, then any other), or None."""
    prefix = Path(mono).resolve().parent.parent
    for gac in (prefix / "lib" / "mono" / "gac" / "Mono.Cecil", Path("/usr/lib/mono/gac/Mono.Cecil")):
        if not gac.is_dir():
            continue
        found = sorted(gac.glob("*/Mono.Cecil.dll"), key=lambda p: (not p.parent.name.startswith(CECIL_VERSION + "__"), p.parent.name))
        if found:
            return found[0]
    return None


def apply_thread_fix(install):
    """Mark FragmentPeakMatcher.maxBondDissociationEnergy [ThreadStatic] in ``install``/StructureFinder.dll (see msfinder_engine).

    Does nothing when the fix is already there. Only the release file (pinned SHA-256) is changed; it is kept as
    StructureFinder.dll.original and the new file replaces it in one rename. Returns the state afterwards ("fixed")."""
    install = Path(install)
    dll = install / mf.STRUCTUREFINDER_DLL
    state = mf.thread_fix_state(install)
    if state == "fixed":
        return state
    if state != "release":
        raise SetupError(f"{dll} is {'missing' if state == 'missing' else 'not the MS-FINDER v3.73 release file (unknown SHA-256)'}; "
                         f"the thread fix is only applied to the release file. Run `{mf.SETUP_COMMAND} --force` to reinstall.")
    mono = mf.find_mono()
    csc = find_csc(mono) if mono else None
    cecil = find_cecil(mono) if mono else None
    if not csc or not cecil:
        raise SetupError("The thread fix needs Mono's C# compiler (csc) and Mono.Cecil from Mono's GAC, which were not found; install the full "
                         "Mono package (for example `brew install mono`).")
    build = work_dir(install) / "threadfix"
    build.mkdir(parents=True, exist_ok=True)
    patcher = build / "ThreadStaticPatch.exe"
    proc = subprocess.run([csc, "-nologo", "-target:exe", f"-out:{patcher}", f"-r:{cecil}", str(PATCH_SOURCE)],
                          capture_output=True, text=True, timeout=600)
    if proc.returncode != 0 or not patcher.is_file():
        raise SetupError("csc failed on the thread fix:\n" + ((proc.stdout or "") + (proc.stderr or ""))[-3000:])
    shutil.copy2(cecil, build / "Mono.Cecil.dll")                   # next to the patcher, whatever the GAC lookup does
    new = install / (mf.STRUCTUREFINDER_DLL + ".threadfix.part")
    new.unlink(missing_ok=True)
    proc = subprocess.run([mono, str(patcher), str(dll), str(new)], capture_output=True, text=True, timeout=600)
    if proc.returncode != 0 or not new.is_file():
        new.unlink(missing_ok=True)
        raise SetupError("The thread fix failed:\n" + ((proc.stdout or "") + (proc.stderr or ""))[-3000:])
    fixed_sha = sha256_of(new)
    original = install / (mf.STRUCTUREFINDER_DLL + ".original")
    if not original.is_file() or sha256_of(original) != mf.STRUCTUREFINDER_SHA256:
        shutil.copy2(dll, original)
    (install / mf.THREAD_FIX_MARKER).write_text(json.dumps({
        "fix": "FragmentPeakMatcher.maxBondDissociationEnergy marked [ThreadStatic] (mzdesk/lib/msfinder/ThreadStaticPatch.cs)",
        "release_sha256": mf.STRUCTUREFINDER_SHA256, "fixed_sha256": fixed_sha, "cecil": str(cecil),
        "same_as_reference_build": fixed_sha == mf.STRUCTUREFINDER_FIXED_SHA256}, indent=1) + "\n", encoding="utf-8")
    os.replace(new, dll)
    state = mf.thread_fix_state(install)
    if state != "fixed":
        raise SetupError(f"The thread fix was written but {dll} does not check out as fixed ({state}).")
    note = "" if fixed_sha == mf.STRUCTUREFINDER_FIXED_SHA256 else " (bytes differ from mzDesk's reference build: another Mono.Cecil version)"
    print(f"Thread fix applied to {dll}{note}; the release file is kept as {original.name}.")
    return state


# ---------------------------------------------------------------- status and setup
def _mono_version(mono):
    try:
        out = subprocess.run([mono, "--version"], capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.TimeoutExpired):
        return "unknown version"
    return (out.splitlines() or ["unknown version"])[0].strip()


def check():
    """Print what is in place; returns 0 when Structure prediction can run, else 1. Starts nothing but `mono --version`."""
    s = mf.status()
    print("MS-FINDER status (mzDesk Structure prediction)")
    print(f"  Mono             : {s['mono'] + '  (' + _mono_version(s['mono']) + ')' if s['mono'] else 'NOT FOUND (install Mono, for example: brew install mono)'}")
    print(f"  Install folder   : {s['dir']}  (MZDESK_MSFINDER_DIR changes it)")
    print(f"  Databases        : {s['resources'] if s['resources'] else 'NOT FOUND (Resources/ is missing)'}")
    exe = s["exe"]
    exe_text = f"{exe}  ({exe.stat().st_size:,} bytes)" if exe else f"NOT FOUND ({mf.EXE_NAME})"
    print(f"  Console exe      : {exe_text}")
    print(f"  Parameter file   : {s['ini'] if s['ini'] else 'NOT FOUND (reinstall mzDesk)'}")
    fix = mf.thread_fix_state(s["dir"])
    print("  Thread fix       : " + {"fixed": "applied", "release": "NOT APPLIED (structure scores change from run to run)",
                                    "missing": "StructureFinder.dll NOT FOUND", "unknown": "UNKNOWN StructureFinder.dll (not the v3.73 release)"}[fix])
    ready = all(s[k] for k in ("mono", "exe", "resources", "ini", "thread_fix"))
    print("  Ready            : " + ("yes" if ready else f"no. Run: {mf.SETUP_COMMAND}"))
    return 0 if ready else 1


def setup(install, force=False, rebuild=False):
    print("Installing MS-FINDER v3.73 (program CC BY 4.0, source LGPL-3.0; the citation is in this module's docstring).")
    mono = mf.find_mono()
    if not mono:
        raise SetupError("Mono was not found. Install it first (macOS: `brew install mono`; Debian / Ubuntu: `apt install mono-complete`), "
                         f"then run `{mf.SETUP_COMMAND}` again.")
    print(f"Mono: {mono} ({_mono_version(mono)})")
    if force or not (install / "Resources").is_dir():
        unzip(get_zip(install, force), install)
        print(f"Unzipped to {install}")
    src_dir = work_dir(install) / "console_src"
    fetch_sources(src_dir, force)
    exe = install / mf.EXE_NAME
    if force or rebuild or not exe.is_file():
        print(f"Compiling {exe}")
        build_console(install, src_dir)
    apply_thread_fix(install)
    print(f"Done: {exe}")
    return check()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Install or check MS-FINDER for the Structure prediction page.")
    ap.add_argument("--check", action="store_true", help="only report what is installed (no download)")
    ap.add_argument("--force", action="store_true", help="download, unzip, fetch and compile again even if present")
    ap.add_argument("--rebuild", action="store_true", help="compile the console exe again from the pinned sources")
    ap.add_argument("--dir", help="install folder (default: MZDESK_MSFINDER_DIR, ~/tools/msfinder/MSFINDER_v3.73 if it exists, else ~/.mzdesk/tools/msfinder/MSFINDER_v3.73)")
    args = ap.parse_args(argv)
    if args.dir:
        os.environ["MZDESK_MSFINDER_DIR"] = args.dir
    if args.check:
        return check()
    try:
        return setup(mf.install_dir(), args.force, args.rebuild)
    except SetupError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
