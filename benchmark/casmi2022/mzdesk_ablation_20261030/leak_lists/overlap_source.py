import csv,json,glob,sys
S='/private/tmp/claude-501/-Users-bowen-Documents-reasearch-01-Core-Projects-mzDesk/46692548-9c05-45e0-abca-108eb47d0813/scratchpad/leak'
R='/path/to/mzDesk'
T=list(csv.DictReader(open(R+'/benchmark/casmi2022/casmi2022_ground_truth.csv')))
ik={int(r['challenge_id']):r['inchikey'][:14] for r in T}
pos={int(r['challenge_id']) for r in T if r['adduct_official'].endswith('+')}
L={}
for f in glob.glob(S+'/lists*.json'): L.update(json.load(open(f)))
def cnt(ids): 
    return len(ids), len([i for i in ids if i in pos]), len([i for i in ids if i not in pos])
def seen(name): s=set(L[name]); return {i for i in ik if ik[i] in s}
res={n:seen(n) for n in L}
for n,v in sorted(res.items()): print(n,len(L[n]),'->',cnt(v))
json.dump({n:sorted(v) for n,v in res.items()},open(S+'/overlap.json','w'))
