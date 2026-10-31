import sys,json
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.inchi import InchiToInchiKey
out={}
for n in ('pos','neg'):
    s=set();bad=0
    for l in open(f'{sys.argv[1]}/{n}.txt'):
        l=l.strip()
        if not l: continue
        k=InchiToInchiKey(l)
        if k: s.add(k[:14])
        else: bad+=1
    out[n]=sorted(s); print(n,len(s),'bad',bad)
json.dump(out,open(sys.argv[2],'w'))
