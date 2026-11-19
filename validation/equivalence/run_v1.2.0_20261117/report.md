# Equivalence experiment: mzDesk route vs direct route

mzDesk 1.2.0 (commit e6937d8), Force merge on, Run FRIGID on, FRIGID on the CPU. Cases: 100 of 100 finished; 0 missing or failed; 89 with every call and every layer candidate identical.

| | Tool level, identical | Layer level, candidates found |
|---|---|---|
| BUDDY | 296 of 296 calls | 632 of 632 |
| Library search | 100 of 100 calls | 123 of 136 (spectra: 670 of 673) |
| MS-FINDER | 93 of 93 runs without a time limit (100 of 106 in all) | 840 of 848 |
| FRIGID | 79 of 79 runs | 185 of 185 |

Library search: compounds by the first 14 InChIKey characters, best score per compound within 0.1; 95 of 100 calls also within 1e-6 (largest best-score difference 0.05002635717391968); the direct search kept 201 spectra that mzDesk did not return; 0 calls had more than mzDesk's 200 returned spectra.
FRIGID devices: cpu.

Operational burden per spectrum (direct route; this script's counting rules): tools 3: 51 spectra, 4: 49 spectra; tool runs 6 (4-9); files written 7 (2-14); parameters set 63 (30-110); hand-offs 6 (1-12) (median, range).
Wall time per case (both routes): 106.5 (34-596) s.

## Calls that differ

- case 028 MS-FINDER msfinder/00 (time-limited): structure files differ: only mzDesk [], only direct ['f0/C30H50O']
- case 029 MS-FINDER msfinder/00 (time-limited): structure files differ: only mzDesk [], only direct ['f0/C25H41NO7', 'f1/C25H38O7']
- case 034 MS-FINDER msfinder/00 (time-limited): f0/C20H32O4: 100 vs 100 structures; best mzDesk [('VNYSSYRCGWBHLG', 7.4401), ('JNUUNUQHXIOFDA', 7.4156), ('BFWYTORDSFIVKP', 7.2548)] vs direct [('VNYSSYRCGWBHLG', 7.4401), ('JNUUNUQHXIOFDA', 7.4156), ('BFWYTORDSFIVKP', 7.2548)]
- case 055 MS-FINDER msfinder/00 (time-limited): f1/C45H82NO8P: 44 vs 45 structures; best mzDesk [('HNFSZKRUPTYBOQ', 6.9752), ('MIPOZDXFSKBSGM', 6.887), ('QCPNXDHEZUMTLB', 6.8821)] vs direct [('HNFSZKRUPTYBOQ', 6.9752), ('MIPOZDXFSKBSGM', 6.887), ('QCPNXDHEZUMTLB', 6.8821)]
- case 063 MS-FINDER msfinder/00 (time-limited): structure files differ: only mzDesk ['f0/C41H78NO9P'], only direct []
- case 080 MS-FINDER msfinder/00 (time-limited): f0/C20H32O3: 90 vs 83 structures; best mzDesk [('KGIJOOYOSFUGPC', 8.2071), ('JSFATNQSLKRBCI', 8.1392), ('DBWQSCSXHFNTMO', 8.0372)] vs direct [('KGIJOOYOSFUGPC', 8.2071), ('JSFATNQSLKRBCI', 8.1392), ('DBWQSCSXHFNTMO', 8.0372)]

## Layer candidates not found

- case 012 Library search: key OVYWMEWYEJLIER, name 6-Quinolinol, accession MSBNK-UFZ-UA006901, library MassBank, score 0.7818882465362549, direct_best 0.7819173336029053
- case 015 Library search: key OYTCZOJKXCTBHG, name Terbutylazine-2-hydroxy, accession MSBNK-Athens_Univ-AU208101, library MassBank, score 0.7356382012367249, direct_best 0.7357413172721863
- case 034 MS-FINDER: ik14 OCBYZNFGCPOXBL, formula C20H32O4, score 7.1679
- case 055 Library search: key GNFCRQUJTIOSKP, name PC(16:0e/15-HEPE), accession MSBNK-RIKEN_IMS-LQA00059, library MassBank, score 0.8571553230285645, direct_best 0.8571568727493286
- case 055 Library search: key ZZJKMJVUYDVCML, name PC(16:0e/11,12-EpETE), accession MSBNK-RIKEN_IMS-LQA00052, library MassBank, score 0.8560906648635864, direct_best 0.8561700582504272
- case 055 Library search: key CQJMWSXXCWAWCC, name PC(16:0e/12-HEPE), accession MSBNK-RIKEN_IMS-LQA00053, library MassBank, score 0.8545243144035339, direct_best 0.8546419739723206
- case 055 Library search: key XHSGENKPTPNVDI, name PC(16:0e/14,15-EpETE), accession MSBNK-RIKEN_IMS-LQA00057, library MassBank, score 0.8298897743225098, direct_best 0.8299345970153809
- case 063 MS-FINDER: ik14 WDDLFSKLSXYHHM, formula C41H78NO9P, score 5.7026
- case 063 MS-FINDER: ik14 BUUGGEDZLZQROW, formula C41H78NO9P, score 5.7011
- case 063 MS-FINDER: ik14 FQCCOAVFQZPNLG, formula C41H78NO9P, score 5.7011
- case 063 MS-FINDER: ik14 QTCAOKZUEKULAU, formula C41H78NO9P, score 5.7006
- case 063 MS-FINDER: ik14 SRFZSYZQMWQOSN, formula C41H78NO9P, score 5.6991
- case 063 MS-FINDER: ik14 UKNPYJPNCNPRBE, formula C41H78NO9P, score 5.6991
- case 066 Library search: key UPFGBRNWZOZZLF, name PC(16:0e/8,9-EET), accession MSBNK-RIKEN_IMS-LQA00071, library MassBank, score 0.8672707080841064, direct_best 0.8672817349433899
- case 066 Library search: key NMKCCEHYPWXWLG, name PC(16:0e/14,15-EET), accession MSBNK-RIKEN_IMS-LQA00056, library MassBank, score 0.7649375796318054, direct_best 0.7649428248405457
- case 066 Library search: key QZDFJRZZGKPADY, name PC(16:0e/15-HETE), accession MSBNK-RIKEN_IMS-LQA00060, library MassBank, score 0.7334724068641663, direct_best 0.7834987640380859
- case 066 Library search: key BXWZXPGVAXFAJU, name PC(16:0e/20-HETE), accession MSBNK-RIKEN_IMS-LQA00067, library MassBank, score 0.7069441080093384, direct_best 0.707267165184021
- case 079 Library search: key WJJMNDUMQPNECX, name Massbank:LU090153 2,6-Pyridinedicarboxylic acid M-H, accession CCMSLIB00005766680, library GNPS, score 0.9221497178077698, direct_best 0.9987831115722656
- case 079 Library search: key GJAWHXHKYYXBSV, name Massbank:LU093053 2,3-Pyridinedicarboxylic acid|Quinolinic acid M-H, accession CCMSLIB00005760333, library GNPS, score 0.9221497178077698, direct_best 0.9987831115722656
- case 080 MS-FINDER: ik14 FVXUUVSFGPLNBJ, formula C20H32O3, score 7.5327
- case 082 Library search: key XWBWIAOWSABHFI, name Cer[NS] d38:1, accession MSBNK-RIKEN_IMS-LQB00425, library MassBank, score 0.8654600381851196, direct_best 0.8657946586608887

