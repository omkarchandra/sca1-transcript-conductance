"""Core representation from frozen S06 transforms; no genotype-coloured outputs."""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

from pathlib import Path
from datetime import datetime,timezone
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator,eigsh
from s06_build import ROOT,A,sha,seal,verify


def weighted_pca(matrix,mean,scale,weights,max_pcs=30):
    scaled=matrix.multiply(1/scale).tocsr();center=mean/scale
    def matvec(v):
        return np.asarray(scaled.T@(weights*(scaled@v))).ravel()-center*np.dot(center,v)
    n=scaled.shape[1]
    operator=LinearOperator((n,n),matvec=matvec,dtype=np.float64)
    if n<=max_pcs+1:
        covariance=np.column_stack([matvec(np.eye(n)[:,i]) for i in range(n)])
        values,vectors=np.linalg.eigh(covariance)
    else:values,vectors=eigsh(operator,k=max_pcs,which='LA',tol=1e-8,v0=np.random.default_rng(20260930).normal(size=n))
    order=np.argsort(values)[::-1];values=values[order];vectors=vectors[:,order]
    retained=(values>max(values[0],0)*1e-12)&(values>0)
    values=values[retained][:max_pcs];vectors=vectors[:,retained][:,:max_pcs]
    for j in range(vectors.shape[1]):
        if vectors[np.argmax(np.abs(vectors[:,j])),j]<0:vectors[:,j]*=-1
    return values,vectors


def main():
    verify(ROOT/_release_paths.expand('${DATA_ROOT}/results/s07_full_contract_frozen_v01'))
    check=ROOT/_release_paths.expand('${DATA_ROOT}/results/s07_workflow_integration_checks_v01/CHECKS.json')
    assert json.loads(check.read_text())['passed']
    verify(A)
    out=ROOT/_release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01');out.mkdir(exist_ok=False)
    normalized=sp.load_npz(A/'normalized_log10k_features.npz')
    cells=pd.read_csv(A/'cell_metadata.tsv.gz',sep='\t');units=pd.read_csv(A/'unit_metadata.tsv',sep='\t',dtype={'animal_id':str})
    keys=json.loads((A/'transform_keys.json').read_text());parameters=np.load(A/'transform_parameters.npz')
    unit=cells.unit_index.to_numpy();records=[]
    for k,key in enumerate(keys):
        mode=key['mode'];pop=key['population'];p_mask=units.cell_type_raw.eq(pop).to_numpy()&(units.n_nuclei>=20).to_numpy()
        fitting=p_mask.copy()
        if mode!='discovery':fitting&=units.fold.ne(int(mode[-1])).to_numpy()
        by_age=units[fitting].groupby('age_weeks').size()
        if len(by_age)<4 or by_age.min()<2:
            records.append(dict(**key,status='insufficient_fitting_mice_per_age'));continue
        eligible=parameters['eligible'][k];mean=parameters['mean'][k,eligible];scale=parameters['scale'][k,eligible]
        fit_cells=np.flatnonzero(fitting[unit]);all_cells=np.flatnonzero(p_mask[unit])
        nfit=int(fitting.sum());weights=1/(units.n_nuclei.to_numpy()[unit[fit_cells]]*nfit)
        matrix=normalized[fit_cells][:,eligible]
        assert np.isclose(weights.sum(),1)
        empirical=np.asarray(matrix.T@weights).ravel();assert np.allclose(empirical,mean,atol=1e-10)
        values,rotation=weighted_pca(matrix,mean,scale,weights)
        if not len(values):records.append(dict(**key,status='no_positive_training_rank'));continue
        pcs=np.asarray(normalized[all_cells][:,eligible].multiply(1/scale)@rotation)-(mean/scale)@rotation
        target=out/(mode+'__'+pop.replace(' ','_'));target.mkdir()
        np.savez_compressed(target/'pca_model.npz',feature_indices=np.flatnonzero(eligible),mean=mean,scale=scale,eigenvalues=values,rotation=rotation)
        np.save(target/'nucleus_PCs.npy',pcs)
        cells.iloc[all_cells].to_csv(target/'nucleus_metadata.tsv.gz',sep='\t',index=False,compression='gzip')
        m=min(200,int(units.loc[fitting,'n_nuclei'].min()));local_lookup=np.full(len(cells),-1,dtype=int);local_lookup[all_cells]=np.arange(len(all_cells))
        landmark=[]
        # Dedicated fixed seed per transform: selection doesn't depend on skips,
        # gene panels or model success. The same IDs later serve broad baseline.
        sampling=np.random.default_rng(20260930)
        for u in np.flatnonzero(fitting):
            selected=sampling.choice(np.flatnonzero(unit==u),size=m,replace=False);landmark.extend(selected.tolist())
        landmark=np.asarray(landmark,dtype=int);lm=cells.iloc[landmark].copy().reset_index(drop=True)
        for field in ['animal_id','age_weeks','genotype','sex']:lm[field]=units[field].to_numpy()[lm.unit_index]
        lm.to_csv(target/'landmark_metadata.tsv',sep='\t',index=False)
        frame=pd.DataFrame(pcs[local_lookup[landmark]],index=lm.obs_name,columns=[f'PC{i+1}' for i in range(pcs.shape[1])]);frame.index.name='obs_name'
        frame.to_csv(target/'landmark_PCs.tsv',sep='\t')
        records.append(dict(**key,status='prepared',n_features=int(eligible.sum()),n_PCs=len(values),n_fitting_mice=nfit,
                            landmarks_per_mouse=m,n_landmarks=len(landmark),no_genotype_coloured_output_viewed=True))
        print('Core PCA prepared:',mode,pop,flush=True)
    (out/'REPRESENTATIONS.json').write_text(json.dumps(records,indent=2)+'\n')
    seal(out,dict(completed_utc=datetime.now(timezone.utc).isoformat(),representations=len(records),
        training_only=True,source_transform_registry_sha256=sha(A/'SHA256SUMS.txt'),
        scientific_tests_fit=False,broad_expression_baseline='same rules/landmark IDs; pending full-transcriptome representation',
        state_registry_gate='no genotype-coloured plots or occupancy before hashed complete state registry',
        code_sha256=sha(Path(__file__)),integration_check_sha256=sha(check)))
    (out/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(out).as_posix()}\n' for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'))


if __name__=='__main__':main()

#%% debug
