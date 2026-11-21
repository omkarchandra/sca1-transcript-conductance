"""Phase 4 Track 2 : the gene -> model-conductance table of the specification, generated from the frozen module definitions (definitions only, no expression data).

Inputs (pinned): ${DATA_ROOT}/s08_module_ids/MODULE_IDS.json (sha256 458713e72b81...) for the stable IDs of modules M1-M5, and ${DATA_ROOT}/results/s10_inputs_v01/GENE_IDENTITIES.tsv
(sha256 e2f3adb9...) for the dataset symbols. The family assignments below are the PROPOSAL of the specification, written as data; every symbol must belong to the stated module (asserted).
Output: one new markdown file (exclusive-create) with the mapping table, the module genes that map to no conductance, and the auxiliary-subunit list. Prints its sha256.
Usage: python ${REPO_ROOT}/05_model_forrest/p4_track2_gene_tables.py <new output file .md>
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

from pathlib import Path
import csv
import hashlib
import json
import sys

ROOT=_ReleasePath(_release_paths.expand('${REPO_ROOT}'))
MODULE_IDS=ROOT/_release_paths.expand('${DATA_ROOT}/s08_module_ids/MODULE_IDS.json')
GENES=ROOT/_release_paths.expand('${DATA_ROOT}/results/s10_inputs_v01/GENE_IDENTITIES.tsv')
PINS={MODULE_IDS:'458713e72b81097063d8b2d51bd784b11b125acee614e65388a815dee15b8866',GENES:'e2f3adb9b5ef0a66e0e002284bbcf1734bdd31816977f661348abb7e380e5a2a'}

# (conductance id, reference-model mechanism(s) in ModelDB 180789 and what they are, module, pore-forming/catalytic genes, auxiliary genes reported but not scaled, certainty note)
FAMILIES=[
 ('g_NaRsg','narsg (soma): resurgent Na+ current','M2',['Scn1a','Scn2a','Scn3a','Scn8a'],['Scn1b','Scn2b','Scn3b','Scn4b'],'mechanism identity from its name; Scn4b is the blocking-particle candidate for the resurgent component (Grieco 2005), reported, not scaled'),
 ('g_h','hpkj (soma), cha (dendrite): HCN current','M2',['Hcn1','Hcn2','Hcn3','Hcn4'],[],'mechanism identity from its name'),
 ('g_BK','bkpkj (soma), kc3 (dendrite): BK-type Ca2+-activated K+','M1',['Kcnma1'],[],'identity from the BK name of bkpkj; kc3 identity UNVERIFIED'),
 ('g_SK','gkca (soma), k23 (dendrite): small-conductance Ca2+-activated K+','M1',['Kcnn1','Kcnn2','Kcnn3'],[],'identity of gkca and k23 UNVERIFIED (Kcnn4 / KCa3.1 is not in M1: not mapped)'),
 ('g_Kv3','kpkj (soma): fast delayed-rectifier-type K+','M1',['Kcnc1','Kcnc2','Kcnc3','Kcnc4'],[],'MODELLING DECISION: gene family of the three Purkinje K+ mechanisms kpkj, kpkj2, kpkjslow is UNVERIFIED; Kv3 for kpkj is a proposal'),
 ('g_Kv1','kv1 (dendrite): Kv1.2-type K+','M1',['Kcna1','Kcna2','Kcna4','Kcna6'],[],'Kcna5 (Kv1.5) is not in M1'),
 ('g_KA','ka (dendrite): A-type K+','M1',['Kcnd1','Kcnd2','Kcnd3'],[],'identity from the name ka'),
 ('g_KM','km (dendrite): M-type K+','M1',['Kcnq2','Kcnq3','Kcnq5'],[],'Kcnq1 and Kcnq4 left unmapped (not M-channel subunits in this proposal)'),
 ('g_KD','kd, khh (dendrite): delayed-rectifier K+','M1',['Kcnb1','Kcnb2'],[],'UNVERIFIED family (proposal: Kv2 for both kd and khh); a gene is assigned to one conductance only'),
 ('g_ERG','erg (dendrite, gbar = 0 in the reference)','M1',['Kcnh2'],[],'held at 0 as in the reference model: the ratio has no effect'),
 ('g_CaP','captain (soma), cap2 (dendrite): P-type Ca2+','M3',['Cacna1a'],[],'identity from the P-type name'),
 ('g_CaT','cat (dendrite): T-type Ca2+','M3',['Cacna1g','Cacna1h','Cacna1i'],[],'identity from the T-type name'),
 ('g_CaR','cae (dendrite): class-E Ca2+','M3',['Cacna1e'],[],'identity from the class-E name'),
 ('p_NaK','pump, mypump, mypumpsoma, newnew: Na+/K+ pump (3 Na+ out : 2 K+ in)','M5',['Atp1a1','Atp1a2','Atp1a3','Atp1a4'],['Atp1b1','Atp1b2','Atp1b3','Fxyd1','Fxyd4','Fxyd5','Fxyd6','Fxyd7'],'alpha = catalytic subunit scaled; beta and FXYD are auxiliary: reported, sensitivity only'),
 ('p_NCX','myexchanger, myexchangersoma: Na+/Ca2+ exchanger','M4',['Slc8a1','Slc8a2','Slc8a3'],[],'Slc8b1 (NCLX) is a mitochondrial exchanger and is not mapped'),
 ('b_Ca','cad, cadiff: Ca2+ buffering and diffusion','M4',['Calb1','Pvalb'],[],'sensitivity analysis only, not in the primary scaling'),
]
HELD_FIXED=[('g_Kslow','kpkj2, kpkjslow (soma): further delayed-rectifier-type K+','the gene family of these two Purkinje K+ mechanisms is UNVERIFIED, so no gene is assigned in the primary analysis; sensitivity: lump with g_Kv3'),('leak','lkpkj (soma), leak (dendrite)','candidate genes in M1/M2 (Kcnk*, Nalcn) are not scaled in the primary analysis'),('kdyn','[K+]o dynamics','no gene'),('nadifl','[Na+]i dynamics','no gene')]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main(argv=None):
    out=Path((argv or sys.argv[1:])[0])
    for path,pin in PINS.items():
        if sha(path)!=pin:raise ValueError(f'input differs from its pin: {path.name}')
    modules=json.loads(MODULE_IDS.read_text())['modules']
    symbol={r['stable_id']:r['dataset_symbol'] for r in csv.DictReader(GENES.open(),delimiter='\t')}
    by_symbol={}
    for module in ('M1','M2','M3','M4','M5'):
        for stable_id in modules[module]['stable_ids']:by_symbol.setdefault(symbol[stable_id],[]).append((module,stable_id))
    for symbol_name,entries in by_symbol.items():
        if len(entries)!=1:raise ValueError(f'symbol in several modules: {symbol_name}')
    def ident(name,module):
        entry=by_symbol.get(name)
        if not entry or entry[0][0]!=module:raise ValueError(f'{name} is not in {module}')
        return f'{name} ({entry[0][1]})'
    lines=['| conductance | reference-model mechanism(s) | module | pore-forming / catalytic genes (scaled) | auxiliary genes (reported, not scaled) | certainty |','|---|---|---|---|---|---|']
    used=set()
    for cid,mech,module,pore,aux,note in FAMILIES:
        used.update(pore)
        lines.append(f"| {cid} | {mech} | {module} | {', '.join(ident(g,module) for g in pore)} | {', '.join(ident(g,module) for g in aux) or '-'} | {note} |")
    aux_all={g for *_,aux,_ in FAMILIES for g in aux}
    unmapped={}
    for module in ('M1','M2','M3','M4','M5'):
        names=sorted(symbol[s] for s in modules[module]['stable_ids'])
        unmapped[module]=[g for g in names if g not in used and g not in aux_all]
    text=[_release_paths.expand('## Gene to conductance table (generated by ${REPO_ROOT}/05_model_forrest/p4_track2_gene_tables.py from MODULE_IDS 458713e7 and GENE_IDENTITIES e2f3adb9)'),'']+lines+['',
        '### Reference-model conductances held fixed','','| item | mechanism | note |','|---|---|---|']+[f'| {a} | {b} | {c} |' for a,b,c in HELD_FIXED]+['','### Module genes that map to no model conductance (not scaled; listed, never substituted)','']
    for module in ('M1','M2','M3','M4','M5'):
        text.append(f"- {module} ({len(modules[module]['stable_ids'])} resolved IDs): "+(', '.join(ident(g,module) for g in unmapped[module]) or 'none'))
    text+=['','### Model channels with no gene in M1-M5 (held fixed)','','- Kv1.5 (Kcna5) and KCa3.1 (Kcnn4), which the Masoli 2015 model contains, have no gene in M1; the Forrest surrogate has no mechanism that needs them.','']
    shared=sorted({g for a in FAMILIES for b in FAMILIES if a[0]<b[0] for g in set(a[3])&set(b[3])})
    text+=['### Genes assigned to more than one conductance (a decision is pending before any run)','',', '.join(shared) or 'none','']
    with out.open('x') as stream:stream.write('\n'.join(text))
    print(hashlib.sha256(out.read_bytes()).hexdigest())

if __name__=='__main__':main()

#%% debug
