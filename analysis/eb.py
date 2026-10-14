"""Empirical Bayes deconvolution of effect sizes: theta/sigma ~ pi0*delta0 + pi1*N(0,t1^2) + pi2*N(0,t2^2), per order."""
import numpy as np, pandas as pd, pickle, warnings; warnings.filterwarnings('ignore')
import datasets as dsm
from lib import truth_fit, Dataset
def fit_mix(z, s, iters=500):
    # z: observed coef/sigma, s: its sd (se/sigma). Components: 0, N(0,t1), N(0,t2)
    p=np.array([.8,.15,.05]); t=np.array([0., .3, 1.5])
    for _ in range(iters):
        v=s[:,None]**2+t[None,:]**2
        L=p[None,:]*np.exp(-z[:,None]**2/(2*v))/np.sqrt(2*np.pi*v); R=L/L.sum(1,keepdims=True)
        p=R.mean(0)
        for k in (1,2):  # moment-style update for slab variance
            w=R[:,k]; t[k]=np.sqrt(max(((w*(z**2-s**2)).sum()/max(w.sum(),1e-9)),1e-4))
        if t[1]>t[2]: t[[1,2]]=t[[2,1]]; p[[1,2]]=p[[2,1]]
    ll=np.log((p[None,:]*np.exp(-z[:,None]**2/(2*(s[:,None]**2+t**2)))/np.sqrt(2*np.pi*(s[:,None]**2+t**2))).sum(1)).sum()
    return p,t,ll
out=[]; fits={}
def run(name, ds, tr, mask=None):
    n=len(ds.cell_of_obs)
    order=np.array([len(t) for t in tr['terms']])
    Z=(tr['coef']/tr['sigma'][None,:])       # T x G, in units of replicate SD (+/-1 coding)
    S=(tr['se']/tr['sigma'][None,:])
    cols=np.arange(Z.shape[1]) if mask is None else mask
    for o in sorted(set(order)):
        lab = {1:'simples',2:'paires',3:'triplets'}.get(o,'ordre 4-5')
        sel=np.isin(order,[o] if o<=3 else [4,5])
        if o==5: continue
        z=Z[sel][:,cols].ravel(); s=S[sel][:,cols].ravel()
        p,t,ll=fit_mix(z,s)
        fits[(name,lab)]=(p,t)
        # in "difference present-absent" units = 2*theta
        out.append(dict(jeu=name,ordre=lab,n_coef=len(z),p_nul=p[0],p_petit=p[1],taille_petit=2*t[1],p_grand=p[2],taille_grand=2*t[2],
                        p_non_nul=p[1]+p[2], frac_grand_parmi_non_nul=p[2]/(p[1]+p[2])))
ds=dsm.ristova(); tr=truth_fit(ds); run('Ristova',ds,tr)
ds=dsm.kellermeier_tx(); tr=truth_fit(ds); run('Kellermeier gènes',ds,tr)
ds=dsm.kellermeier(); tr=truth_fit(ds); run('Kellermeier traits',ds,tr)
ds=dsm.tekin(); sub=np.sort(np.random.default_rng(7).choice(ds.Y.shape[1],2000,replace=False))
ds=Dataset(ds.name,ds.factors,ds.cells,ds.cell_of_obs,ds.Y[:,sub],[ds.resp_names[i] for i in sub]); tr=truth_fit(ds); run('Tekin',ds,tr)
df=pd.DataFrame(out); pd.set_option('display.width',200); print(df.round(3).to_string(index=False))
df.to_csv('fig/effect_distributions.csv',index=False); pickle.dump(fits,open('prep/eb_fits.pkl','wb'))
