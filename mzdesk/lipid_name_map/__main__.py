"""``python3 -m lipid_name_map``：命令行入口，用法见 ``lipid_name_map.cli`` 与项目根的 README.md。"""
from __future__ import annotations

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
