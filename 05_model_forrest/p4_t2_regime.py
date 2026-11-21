"""Phase 4 Track 2 : regime classifier for a settled somatic voltage window (spec TRACK2_SPEC_v01 sections 4.1 and 4.2; v02 of the module: contiguous segmentation and one edge rule). Pure numpy; knows nothing about NEURON.

Labels: TONIC, BURSTING, BIMODAL (two of tonic / bursting / quiescence alternate), TRIMODAL (all three), QUIESCENT, DEPOLARIZATION_BLOCK, OTHER.
BIMODAL is added to the spec's list so that the bimodal pattern of the acceptance gate has its own label; the spec's label list is amended in v02. A window that does not meet a definition is OTHER, never forced into a class.
All numeric bounds are arguments (class Bounds), to be fixed ONCE from the reference runs before any mouse is simulated.
Definitions (window of at least min_window_ms after settling):
- spikes: upward crossings of spike_threshold_mV, at least refractory_ms apart.
- no spikes: QUIESCENT if mean V < v_quiescent_max and the mean V of the two halves differs by at most drift_max; DEPOLARIZATION_BLOCK if mean V > v_block_min and the drift is at most drift_max; otherwise OTHER.  1 or 2 spikes: OTHER.
- gaps (quiescent episodes): every interval between spikes above gap_silent_ms, and the leading silence (window start to first spike) and trailing silence (last spike to window end) when they are at least gap_silent_ms. One edge rule: a silence at a
  window edge counts like any other gap, and an epoch is 'cut by the window' (touches the edge) only when it starts or ends within gap_silent_ms of the window edge, so a gap is never counted on one side of an ignored fragment and ignored on the other.
- epochs: the spike train split at the internal gaps. An epoch with fewer than min_epoch_spikes spikes is SHORT. Longer epochs are segmented WITHOUT needing a silence: the CV of every run of block_isis consecutive intervals is labelled TONIC (below cv_tonic_max),
  BURSTING (at least cv_burst_min) or AMBIGUOUS, clean runs shorter than block_isis//2 windows are relabelled AMBIGUOUS, and an AMBIGUOUS run is a transition only if it is shorter than block_isis windows and its two neighbours differ (two clean runs of different kinds, or a clean run on one side and the epoch end on the other; an AMBIGUOUS run with the same kind on both sides, or the whole epoch, is not a transition); any other AMBIGUOUS run, and any epoch shorter than block_isis intervals whose CV is in between, is an AMBIGUOUS segment.
- an interior SHORT or AMBIGUOUS segment makes the window OTHER; SHORT epochs cut by a window edge are ignored (their adjacent gap is still counted). Modes = kinds of the clean segments plus QUIESCENT when a gap exists;
  one mode with no gap gives TONIC or BURSTING, two modes BIMODAL, three TRIMODAL.
Usage: from p4_t2_regime import Bounds, classify
"""
from dataclasses import dataclass
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

LABELS=('TONIC','BURSTING','BIMODAL','TRIMODAL','QUIESCENT','DEPOLARIZATION_BLOCK','OTHER')

@dataclass(frozen=True)
class Bounds:
    v_quiescent_max:float
    v_block_min:float
    spike_threshold_mV:float=0.0
    refractory_ms:float=1.0
    drift_max:float=0.1
    gap_silent_ms:float=1000.0
    min_epoch_spikes:int=5
    cv_tonic_max:float=0.3
    cv_burst_min:float=0.5
    block_isis:int=8
    min_window_ms:float=20000.0

    def __post_init__(self):
        if self.block_isis<4:raise ValueError('block_isis must be at least 4 (the smoothing length block_isis//2 must be at least 2)')

def detect_spikes(t,v,bounds):
    t=np.asarray(t,float);v=np.asarray(v,float)
    up=np.flatnonzero((v[:-1]<bounds.spike_threshold_mV)&(v[1:]>=bounds.spike_threshold_mV))
    times=t[up]+(bounds.spike_threshold_mV-v[up])/(v[up+1]-v[up])*(t[up+1]-t[up])
    kept=[]
    for time in times:
        if not kept or time-kept[-1]>=bounds.refractory_ms:kept.append(time)
    return np.array(kept)

def label_cv(cv,bounds):
    if cv<bounds.cv_tonic_max:return 'TONIC'
    if cv>=bounds.cv_burst_min:return 'BURSTING'
    return 'AMBIGUOUS'

def segment_epoch(times,bounds):
    """Clean kinds and ambiguity of one epoch (spike times), in order: a list of kinds from TONIC, BURSTING, AMBIGUOUS."""
    isi=np.diff(times);k=bounds.block_isis
    if len(isi)<k:return [label_cv(float(isi.std()/isi.mean()),bounds)]
    windows=sliding_window_view(isi,k)
    labels=[label_cv(float(w.std()/w.mean()),bounds) for w in windows]
    runs=[]
    for label in labels:
        if runs and runs[-1][0]==label:runs[-1][1]+=1
        else:runs.append([label,1])
    for run in runs:
        if run[1]<k//2:run[0]='AMBIGUOUS'
    merged=[]
    for run in runs:
        if merged and merged[-1][0]==run[0]:merged[-1][1]+=run[1]
        else:merged.append(list(run))
    out=[]
    for i,(label,length) in enumerate(merged):
        if label!='AMBIGUOUS':out.append(label);continue
        before=merged[i-1][0] if i>0 else None;after=merged[i+1][0] if i+1<len(merged) else None
        transition=length<k and before!=after
        if not transition:out.append('AMBIGUOUS')
    return out

def epochs(spikes,t0,t1,bounds):
    intervals=np.diff(spikes)
    cut=np.flatnonzero(intervals>bounds.gap_silent_ms)
    groups=np.split(np.arange(len(spikes)),cut+1)
    gaps=len(cut)+int(spikes[0]-t0>=bounds.gap_silent_ms)+int(t1-spikes[-1]>=bounds.gap_silent_ms)
    out=[]
    for k,group in enumerate(groups):
        times=spikes[group];n=len(times)
        edge=(k==0 and times[0]-t0<bounds.gap_silent_ms) or (k==len(groups)-1 and t1-times[-1]<bounds.gap_silent_ms)
        kinds=['SHORT'] if n<bounds.min_epoch_spikes else segment_epoch(times,bounds)
        out.append({'start':float(times[0]),'end':float(times[-1]),'n':int(n),'kinds':kinds,'touches_window_edge':bool(edge)})
    return out,gaps

def classify(t,v,bounds):
    t=np.asarray(t,float);v=np.asarray(v,float)
    span=float(t[-1]-t[0])
    info={'window_ms':span,'mean_v':float(v.mean())}
    if span<bounds.min_window_ms:return dict(info,label='OTHER',reason='window shorter than min_window_ms')
    half=len(v)//2
    drift=abs(float(v[half:].mean()-v[:half].mean()))
    spikes=detect_spikes(t,v,bounds)
    info.update(drift=drift,n_spikes=int(len(spikes)),rate_hz=float(len(spikes)/span*1000))
    if len(spikes)==0:
        if drift>bounds.drift_max:return dict(info,label='OTHER',reason='no spikes, V not settled')
        if info['mean_v']<bounds.v_quiescent_max:return dict(info,label='QUIESCENT',reason='no spikes, low V')
        if info['mean_v']>bounds.v_block_min:return dict(info,label='DEPOLARIZATION_BLOCK',reason='no spikes, high V')
        return dict(info,label='OTHER',reason='no spikes, V between the bounds')
    if len(spikes)<3:return dict(info,label='OTHER',reason='fewer than 3 spikes')
    parts,gaps=epochs(spikes,float(t[0]),float(t[-1]),bounds)
    info.update(epochs=parts,n_gaps=gaps)
    kept=[]
    for part in parts:
        if part['kinds']==['SHORT'] and part['touches_window_edge']:continue
        for kind in part['kinds']:
            if kind in ('SHORT','AMBIGUOUS'):return dict(info,label='OTHER',reason=f'interior {kind} segment')
            kept.append(kind)
    if not kept:return dict(info,label='OTHER',reason='no usable epoch')
    modes=set(kept)|({'QUIESCENT'} if gaps else set())
    info['modes']=sorted(modes)
    if len(modes)==1:return dict(info,label=kept[0],reason='single mode')
    return dict(info,label='BIMODAL' if len(modes)==2 else 'TRIMODAL',reason='several modes')

#%% debug
