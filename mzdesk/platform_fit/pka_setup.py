"""MolGpKa（pKa 预测）权重准备，一次性运行：python -m mzdesk.platform_fit.pka_setup

MolGpKa 的权重就放在作者的 GitHub 仓库里（Xundrug/MolGpKa 的 models/，MIT 协议）。这里按固定的提交号下载
酸、碱两个 .pth（各约 6 MB），校验 sha256，再转成 numpy 的 molgpka_weights.npz。只有这一步需要 PyTorch
（读 .pth），运行时不需要 torch，也不需要 torch_geometric。LICENSE 一并下载，放在权重旁边。
文件都放在 ~/.mzdesk/models/molgpka/（环境变量 MZDESK_MOLGPKA_DIR 可改目录）。权重不入库。
"""
import argparse
import hashlib
import os
import sys
import urllib.request
from pathlib import Path

COMMIT = "4dc835219baad18c2e705274ac54d49597703966"      # Xundrug/MolGpKa main，2024-01-11
RAW = f"https://raw.githubusercontent.com/Xundrug/MolGpKa/{COMMIT}/"
FILES = {                                                # 仓库里的路径 -> sha256
    "models/weight_acid.pth": "e0cb05874937f17a4df6943f37beff1c76d14d760787c57d553ac56390147b39",
    "models/weight_base.pth": "212021f43168208bd67a5c15750f9a176cf4d5906d73a85718f4b15ff8d4b95d",
}
LICENSE_FILE = "LICENSE.md"
LICENSE_SHA256 = None                                    # 许可证文本只存档，不校验
KINDS = ("acid", "base")
SETUP_COMMAND = "python -m mzdesk.platform_fit.pka_setup"


def model_dir():
    return Path(os.environ.get("MZDESK_MOLGPKA_DIR") or Path.home() / ".mzdesk" / "models" / "molgpka")


def default_weights():
    return model_dir() / "molgpka_weights.npz"


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _fetch(rel, dest):
    part = dest.with_suffix(dest.suffix + ".part")
    print(f"Downloading {RAW}{rel}\n  -> {dest}")
    request = urllib.request.Request(RAW + rel, headers={"User-Agent": "mzdesk"})
    with urllib.request.urlopen(request, timeout=60) as resp, open(part, "wb") as out:
        while block := resp.read(1 << 20):
            out.write(block)
    part.rename(dest)


def get_file(rel, sha):
    """Path of the verified file ``rel`` of the pinned commit; downloads it when missing."""
    dest = model_dir() / Path(rel).name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        _fetch(rel, dest)
    if sha:
        digest = sha256_of(dest)
        if digest != sha:
            raise SystemExit(f"{dest} has sha256 {digest}, expected {sha}. Delete it and run again.")
    return dest


def convert(pth_paths, npz_path):
    """{'acid': .pth, 'base': .pth} -> npz with names '<kind>.<state_dict key>' (float32; num_batches_tracked dropped)."""
    import numpy as np

    try:
        import torch
    except ImportError:
        raise SystemExit("PyTorch is needed once to read the original .pth files:\n"
                         "  pip install torch\nThen run this command again. "
                         "(pKa prediction itself does not use torch.)") from None
    tensors = {}
    for kind, path in pth_paths.items():
        state = torch.load(path, map_location="cpu", weights_only=True)      # tensors only: no pickled code is run
        for name, t in state.items():
            if not name.endswith("num_batches_tracked"):
                tensors[f"{kind}.{name}"] = t.numpy().astype("float32")
    tmp = npz_path.with_suffix(".tmp.npz")
    np.savez(tmp, **tensors)
    tmp.replace(npz_path)
    return len(tensors)


def setup(force=False):
    npz = default_weights()
    if npz.exists() and not force:
        print(f"{npz} already exists (use --force to rebuild).")
        return npz
    paths = {kind: get_file(rel, sha) for kind, (rel, sha) in zip(KINDS, FILES.items())}
    get_file(LICENSE_FILE, LICENSE_SHA256)
    n = convert(paths, npz)
    print(f"Wrote {n} tensors to {npz} ({npz.stat().st_size / 1e6:.0f} MB). MolGpKa is MIT-licensed; its LICENSE is next to the weights.")
    return npz


def main(argv=None):
    ap = argparse.ArgumentParser(description="Download the MolGpKa weights and convert them to numpy.")
    ap.add_argument("--force", action="store_true", help="rebuild molgpka_weights.npz even if it exists")
    setup(ap.parse_args(argv).force)


if __name__ == "__main__":
    sys.exit(main())
