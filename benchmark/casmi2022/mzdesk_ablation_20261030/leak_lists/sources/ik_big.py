"""usage: ik_big.py in_file smiles_col(0-based, or -1 = whole line) skip_header out_json name [delimiter]
Only keeps IK14 that are in the 500 truths (memory-light)."""
import sys,csv,json
from multiprocessing import Pool
from rdkit import Chem,RDLogger
from rdkit.Chem.inchi import MolToInchiKey
RDLogger.DisableLog('rdApp.*')
R='/path/to/mzDesk'
TR={r['inchikey'][:14] for r in csv.DictReader(open(R+'/benchmark/casmi2022/casmi2022_ground_truth.csv'))}
def work(lines):
    hit=set();n=0;bad=0
    for s in lines:
        m=Chem.MolFromSmiles(s)
        if m is None: bad+=1;continue
        try: k=MolToInchiKey(m)
        except Exception: bad+=1;continue
        n+=1
        if k and k[:14] in TR: hit.add(k[:14])
    return hit,n,bad
def chunks(f,col,skip,delim,size=20000):
    buf=[]
    for i,l in enumerate(f):
        if i<skip: continue
        l=l.rstrip('\n')
        if not l: continue
        buf.append(l if col<0 else l.split(delim)[col])
        if len(buf)==size: yield buf;buf=[]
    if buf: yield buf
if __name__=='__main__':
    fn,col,skip,out,name=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],sys.argv[5]
    delim=sys.argv[6].replace('\\t','\t') if len(sys.argv)>6 else '\t'
    hit=set();N=0;B=0
    with Pool(8) as p, open(fn) as f:
        for h,n,b in p.imap(work,chunks(f,col,skip,delim)):
            hit|=h;N+=n;B+=b
            if (N//20000)%20==0: print(N,B,len(hit),flush=True)
    json.dump({name:sorted(hit),name+'_n':[N,B]},open(out,'w'));print('done',N,B,len(hit))
