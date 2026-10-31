import csv,json,pickle
from rdkit import Chem,RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.inchi import MolToInchiKey
class U(pickle.Unpickler):
    def find_class(self,m,n): raise pickle.UnpicklingError(m+n)
S='/private/tmp/claude-501/-Users-bowen-Documents-reasearch-01-Core-Projects-mzDesk/46692548-9c05-45e0-abca-108eb47d0813/scratchpad/leak'
R='/path/to/mzDesk'
out={}
# JESTR NPLIB1
d=U(open(R+'/benchmark/JESTR/data/NPLIB1/split.pkl','rb')).load()
out['jestr_nplib1_train']=sorted({k[:14] for k in d['train']})
out['jestr_nplib1_trainval']=sorted({k[:14] for k in d['train']+d['valid']})
# JESTR MassSpecGym (smiles)
d=U(open(R+'/benchmark/JESTR/data/MassSpecGym/split.pkl','rb')).load()
s=set()
for smi in d['train']:
    m=Chem.MolFromSmiles(smi)
    if m: s.add(MolToInchiKey(m)[:14])
out['jestr_msg_train']=sorted(s)
# FRIGID local MSG labels/split: train spectra
lab={r['spec']:r for r in csv.DictReader(open('~/tools/frigid/FRIGID/data/msg/labels.tsv'),delimiter='\t')}
sp=[r for r in csv.DictReader(open('~/tools/frigid/FRIGID/data/msg/split.tsv'),delimiter='\t')]
import collections
print(collections.Counter(r['split'] for r in sp))
out['msg_train']=sorted({lab[r['name']]['inchikey'][:14] for r in sp if r['split']=='train'})
out['msg_all']=sorted({r['inchikey'][:14] for r in lab.values()})
out.update(json.load(open(S+'/ik_sirius.json')))
out={('sirius_'+k if k in('pos','neg') else k):v for k,v in out.items()}
json.dump(out,open(S+'/lists1.json','w'))
print({k:len(v) for k,v in out.items()})
