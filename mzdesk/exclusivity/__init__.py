"""排他性：同分子式候选、同量异位干扰、异构体。A 类。

同分子式与同量异位由 mzdesk.ms1.exact_mass 查 PubChem（在线），在 api/runners/exclusivity.py 里组表；
异构体在 isomers.py（MS-FINDER 本地库）。m/z 关系、代谢物预测、脂质酰基枚举已从页面去掉。
"""
