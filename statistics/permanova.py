import numpy as np, pandas as pd
from scipy.spatial.distance import pdist, squareform
from bio_stats import *
raw = load()
def mat(groups):
    X=[];lab=[];loc=[]
    for g in groups:
        b=bio_matrix(raw,g).T; X.append(b); lab+= [('inf' if g.startswith('infected') else 'non' if g.lower().startswith('nonhost') else 'un')]*len(b); loc+=[g[-3:]]*len(b)
    return np.vstack(X), np.array(lab), np.array(loc)
def pseudoF(D,lab):
    n=len(lab); g=np.unique(lab); A=-0.5*D**2
    SST=(D[np.triu_indices(n,1)]**2).sum()/n
    SSW=sum((D[np.ix_(lab==k,lab==k)][np.triu_indices((lab==k).sum(),1)]**2).sum()/(lab==k).sum() for k in g)
    a=len(g); return ((SST-SSW)/(a-1))/(SSW/(n-a)), (SST-SSW)/SST
def perm(X,lab,strata,metric,nperm=9999,seed=0):
    Xl=np.log10(X+1); keep=Xl.std(0)>0; Xl=Xl[:,keep]
    if metric=='euc': Z=(Xl-Xl.mean(0))/Xl.std(0); D=squareform(pdist(Z))
    else: D=squareform(pdist(X[:,keep],'braycurtis'))
    F0,R2=pseudoF(D,lab); rng=np.random.default_rng(seed); c=0
    for _ in range(nperm):
        l=lab.copy()
        for s in np.unique(strata):
            idx=np.where(strata==s)[0]; l[idx]=rng.permutation(l[idx])
        c+= pseudoF(D,l)[0]>=F0
    return F0,R2,(c+1)/(nperm+1)
tests={'infected vs uninfected (pooled, strata=locality)':(['infectedTHAI','infectedCAM','infectedILO','UNinfectedTHAI','UNinfectedCAM','UNinfectedILO'],True),
       'infected vs uninfected ILO':(['infectedILO','UNinfectedILO'],False),
       'non-host vs infected ILO':(['nonhostILO','infectedILO'],False),
       'non-host vs infected CAM':(['NonhostCAM','infectedCAM'],False),
       'non-host vs infected (CAM+ILO, strata)':(['NonhostCAM','infectedCAM','nonhostILO','infectedILO'],True)}
for k,(gs,st) in tests.items():
    X,lab,loc=mat(gs); strata=loc if st else np.zeros(len(lab))
    for met in ['euc','bc']:
        F,R2,p=perm(X,lab,strata,met, nperm=9999)
        print(f"{k:45s} {met}: n={dict(zip(*np.unique(lab,return_counts=True)))} F={F:.2f} R2={R2:.3f} p={p:.4f}")
