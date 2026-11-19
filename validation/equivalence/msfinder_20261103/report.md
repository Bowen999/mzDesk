# MS-FINDER consistency: mzDesk vs standalone

mzDesk 0.5.6 (commit c32f8f2), MS-FINDER v3.73 with mzDesk's thread fix (~/tools/msfinder/MSFINDER_v3.73).
Cases: 100 of 100 finished; 0 without an MS-FINDER run; 0 missing or failed.

| Measure | Consistency |
|---|---|
| MS-FINDER runs, identical (formulas, structures, scores) | 94/101 = 93.1%  (95% CI 86.4-96.6%) |
| ... without time-limited runs | 87/87 = 100.0%  (95% CI 95.8-100.0%) |
| ... time-limited runs only | 7/14 = 50.0%  (95% CI 26.8-73.2%) |
| MS-FINDER runs, byte-identical structure files | 94/101 = 93.1%  (95% CI 86.4-96.6%) |
| Cases, every run identical | 93/100 = 93.0%  (95% CI 86.3-96.6%) |
| ... without time-limited runs | 86/86 = 100.0%  (95% CI 95.7-100.0%) |
| Cases, mzDesk's MS-FINDER layer fully found in the standalone output | 98/100 = 98.0%  (95% CI 93.0-99.4%) |
| Layer candidates found with the same score | 845/848 = 99.6%  (95% CI 99.0-99.9%) |

## Runs that differ

- case 028 runs/00 (time-limited): f0/C30H50O: 92 vs 92 structures; best mzDesk [('CAHGCLMLTWQZNJ', 8.3394), ('GGGUGZHBAOMSFJ', 8.0565), ('FSLPMRQHCOLESF', 7.9886)] vs direct [('CAHGCLMLTWQZNJ', 8.3394), ('GGGUGZHBAOMSFJ', 8.0565), ('FSLPMRQHCOLESF', 7.9886)]
- case 029 runs/00 (time-limited): f1/C25H38O7: 45 vs 46 structures; best mzDesk [('LYTPVRMVQVQYGM', 5.1583), ('QRWBFJVOYWUREI', 5.1037), ('WWPCIPLSEBENOX', 5.1037)] vs direct [('LYTPVRMVQVQYGM', 5.1583), ('QRWBFJVOYWUREI', 5.1037), ('WWPCIPLSEBENOX', 5.1037)]
- case 034 runs/00 (time-limited): f0/C20H32O4: 100 vs 100 structures; best mzDesk [('VNYSSYRCGWBHLG', 7.4401), ('JNUUNUQHXIOFDA', 7.4156), ('BFWYTORDSFIVKP', 7.2548)] vs direct [('VNYSSYRCGWBHLG', 7.4401), ('JNUUNUQHXIOFDA', 7.4156), ('BFWYTORDSFIVKP', 7.2548)]
- case 055 runs/00 (time-limited): f1/C45H82NO8P: 45 vs 46 structures; best mzDesk [('HNFSZKRUPTYBOQ', 6.9752), ('MIPOZDXFSKBSGM', 6.887), ('QCPNXDHEZUMTLB', 6.8821)] vs direct [('HNFSZKRUPTYBOQ', 6.9752), ('MIPOZDXFSKBSGM', 6.887), ('QCPNXDHEZUMTLB', 6.8821)]
- case 063 runs/00 (time-limited): f0/C41H78NO9P: 63 vs 61 structures; best mzDesk [('BUUGGEDZLZQROW', 5.7011), ('FQCCOAVFQZPNLG', 5.7011), ('QTCAOKZUEKULAU', 5.7006)] vs direct [('WDDLFSKLSXYHHM', 5.7026), ('BUUGGEDZLZQROW', 5.7011), ('FQCCOAVFQZPNLG', 5.7011)]
- case 080 runs/00 (time-limited): f0/C20H32O3: 100 vs 100 structures; best mzDesk [('KGIJOOYOSFUGPC', 8.2071), ('NLUNAYAEIJYXRB', 8.1737), ('JSFATNQSLKRBCI', 8.1392)] vs direct [('KGIJOOYOSFUGPC', 8.2071), ('NLUNAYAEIJYXRB', 8.1737), ('JSFATNQSLKRBCI', 8.1392)]
- case 088 runs/00 (time-limited): f1/C43H78NO8P: 50 vs 48 structures; best mzDesk [('ANRKEHNWXKCXDB', 7.2855), ('SQGZFCFLUVPOSZ', 7.1646), ('MQXKQCUGNLVKIN', 7.1568)] vs direct [('ANRKEHNWXKCXDB', 7.2855), ('SQGZFCFLUVPOSZ', 7.1646), ('MQXKQCUGNLVKIN', 7.1568)]

