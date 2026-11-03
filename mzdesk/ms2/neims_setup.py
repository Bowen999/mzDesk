"""NEIMS（EI 质谱预测）权重准备，一次性运行：python -m mzdesk.ms2.neims_setup

官方权重地址（Google 存储桶）已失效，这里从 Wayback Machine 的存档下载同一个 zip，校验 sha256 后解压，
再把推理需要的张量转成 numpy 的 neims_weights.npz（约 170 MB）。只有这一步需要 TensorFlow，运行时不需要。
文件都放在 ~/.mzdesk/models/neims/（环境变量 MZDESK_NEIMS_DIR 可改目录）。
"""
import argparse
import hashlib
import os
import sys
import urllib.request
import zipfile
from pathlib import Path

ZIP_URL = ("https://web.archive.org/web/20230627165112id_/"
           "https://storage.googleapis.com/deep-molecular-massspec/massspec_weights/massspec_weights.zip")
ZIP_SHA256 = "3c80ee1b7d9738bce5c0a22092728bf17aa1932dcd44859db976ddaa9ed0aa52"
CHECKPOINT = "model.ckpt-100000"
SCOPE = "spectrum_predictor/"         # 只取推理用的变量，跳过 Adam 等优化器状态
SETUP_COMMAND = "python -m mzdesk.ms2.neims_setup"


def model_dir():
    return Path(os.environ.get("MZDESK_NEIMS_DIR") or Path.home() / ".mzdesk" / "models" / "neims")


def default_weights():
    return model_dir() / "neims_weights.npz"


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _download(dest):
    part = dest.with_suffix(".part")
    print(f"Downloading {ZIP_URL}\n  -> {dest} (about 466 MB, this can take a while)")
    request = urllib.request.Request(ZIP_URL, headers={"User-Agent": "mzdesk"})
    with urllib.request.urlopen(request, timeout=60) as resp, open(part, "wb") as out:
        total = int(resp.headers.get("Content-Length") or 0)
        done, report = 0, 10 << 20
        while block := resp.read(1 << 20):
            out.write(block)
            done += len(block)
            if done >= report:
                print(f"  {done >> 20} / {total >> 20} MB" if total else f"  {done >> 20} MB", flush=True)
                report += 10 << 20
    part.replace(dest)


def get_zip():
    """Path of the verified weights zip; downloads it when missing."""
    dest = model_dir() / "massspec_weights.zip"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        _download(dest)
    digest = sha256_of(dest)
    if digest != ZIP_SHA256:
        raise SystemExit(f"{dest} has sha256 {digest}, expected {ZIP_SHA256}. Delete it and run again.")
    return dest


def unzip(zip_path):
    """Extract the checkpoint files (not the .meta graph) next to the zip; returns the checkpoint prefix."""
    out = model_dir()
    want = [f"massspec_weights/{CHECKPOINT}.index", f"massspec_weights/{CHECKPOINT}.data-00000-of-00001"]
    with zipfile.ZipFile(zip_path) as z:
        for name in want:
            if not (out / name).exists():
                z.extract(name, out)
    return out / "massspec_weights" / CHECKPOINT


def convert(prefix, npz_path):
    """Checkpoint -> npz with the original variable names ('/' written as '.'), float32."""
    import numpy as np

    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    try:
        import tensorflow as tf
    except ImportError:
        raise SystemExit("TensorFlow is needed once to read the original checkpoint:\n"
                         "  pip install tensorflow\nThen run this command again. "
                         "(NEIMS prediction itself does not use TensorFlow.)") from None
    reader = tf.train.load_checkpoint(str(prefix))
    names = sorted(n for n in reader.get_variable_to_shape_map() if n.startswith(SCOPE))
    tensors = {n[len(SCOPE):].replace("/", "."): reader.get_tensor(n).astype("float32") for n in names}
    tmp = npz_path.with_suffix(".tmp.npz")
    np.savez(tmp, **tensors)
    tmp.replace(npz_path)
    return len(tensors)


def setup(force=False):
    npz = default_weights()
    if npz.exists() and not force:
        print(f"{npz} already exists (use --force to rebuild).")
        return npz
    prefix = unzip(get_zip())
    n = convert(prefix, npz)
    print(f"Wrote {n} tensors to {npz} ({npz.stat().st_size / 1e6:.0f} MB).")
    return npz


def main(argv=None):
    ap = argparse.ArgumentParser(description="Download the NEIMS weights and convert them to numpy.")
    ap.add_argument("--force", action="store_true", help="rebuild neims_weights.npz even if it exists")
    setup(ap.parse_args(argv).force)


if __name__ == "__main__":
    sys.exit(main())
