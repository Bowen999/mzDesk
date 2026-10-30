"""模块 id -> runner。

runner(inputs: dict) -> {"outputs": {块 id: 数据}, "warnings": [str]}。
spec 里 live=True 的模块都必须在 RUNNERS 里（tests/test_api.py 检查）。引擎在 runner 内部才导入，
所以某个引擎缺依赖时，其他模块照常可用。
"""
from mzdesk.api.runners import (
    adducts, calculator, ccs, ei_ms, esi_msms, exclusivity, gnps, identity, inorganic, ionization, met_shift, mrm, ms1, ms2, ms2_denovo, ms2_structure,
    nca, non_ms, stability, structure,
)

RUNNERS = {**identity.RUNNERS, **structure.RUNNERS, **ms1.RUNNERS, **ms2.RUNNERS, **ms2_structure.RUNNERS, **ms2_denovo.RUNNERS, **ei_ms.RUNNERS, **esi_msms.RUNNERS, **ionization.RUNNERS, **non_ms.RUNNERS,
           **adducts.RUNNERS, **inorganic.RUNNERS, **calculator.RUNNERS, **met_shift.RUNNERS, **mrm.RUNNERS, **nca.RUNNERS, **stability.RUNNERS,
           **exclusivity.RUNNERS, **gnps.RUNNERS, **ccs.RUNNERS}
