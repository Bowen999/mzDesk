"""模块 id -> runner。

runner(inputs: dict) -> {"outputs": {块 id: 数据}, "warnings": [str]}。
spec 里 live=True 的模块都必须在 RUNNERS 里（tests/test_api.py 检查）。引擎在 runner 内部才导入，
所以某个引擎缺依赖时，其他模块照常可用。
"""
from mzdesk.api.runners import identity, ms1, ms2, structure

RUNNERS = {**identity.RUNNERS, **structure.RUNNERS, **ms1.RUNNERS, **ms2.RUNNERS}
