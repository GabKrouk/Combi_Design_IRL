"""Recalibrated planning simulation for 18 compounds: recall of M target interactions of size d (difference
present-absent, in replicate-SD units) among S=14 active compounds, with a background of small interactions."""
import sys, numpy as np, pandas as pd, warnings, itertools; warnings.filterwarnings('ignore')
from multiprocessing import Pool
from lib import Dataset, recover, terms, model_cols
K=18; R=3
def one(args):
    order, d, B, N, seed, G = args
    rng=np.random.default_rng(seed)
    A=np.zeros((N,K),int)
    for j in range(K): A[rng.permutation(N)[:N//2],j]=1
    allt=terms(K,3); o=np.array([len(t) for t in allt]); X=model_cols(np.repeat(A,R,0),allt)
    pairs=np.where(o==2)[0]; trips=np.where(o==3)[0]
    th=np.zeros((len(allt),G))
    for g in range(G):
        m=rng.choice(K,14,replace=False); th[m,g]=rng.normal(0,1.1,14)          # realistic main effects (coef, +/-1 coding)
        tgt=rng.choice(pairs if order==2 else trips,5,replace=False); th[tgt,g]=rng.choice([-1,1],5)*d/2
        if B>0:
            bg=np.setdiff1d(np.concatenate([pairs,trips]),tgt); th[bg,g]+=rng.normal(0,np.sqrt(B/len(bg)),len(bg))
    Y=X@th+rng.normal(size=(X.shape[0],G))
    ds=Dataset('s',list(range(K)),A,np.repeat(np.arange(N),R),Y,range(G))
    res=recover(ds,np.arange(N),cand_order=order,max_steps=15)
    idx2={t:i for i,t in enumerate(terms(K,order))}; full={t:i for i,t in enumerate(allt)}
    rec=[];fdr=[];ok=[]
    for g in range(G):
        tgt={idx2[allt[j]] for j in np.where((np.abs(th[:,g])>=d/2-1e-12)&(o==order))[0]}
        s=set(res['sel'][g]); hit=len(s&tgt); rec.append(hit/5); fdr.append((len(s)-hit)/len(s) if s else 0); ok.append(s==tgt)
    return dict(order=order,d=d,B=B,N=N,recall=np.mean(rec),fdr=np.mean(fdr),exact=np.mean(ok))
if __name__=='__main__':
    jobs=[]; s=0
    for order,Ns in [(2,[40,60,80,120,160,240,320]),(3,[60,80,120,160,240,320,480])]:
        for d in [0.5,1.0,1.5,2.0,3.0]:
            for B in [0.02,0.1,0.45]:
                for N in Ns: s+=1; jobs.append((order,d,B,N,s,60))
    with Pool(2) as p: rows=p.map(one,jobs,chunksize=1)
    df=pd.DataFrame(rows); df.to_csv('fig/sim18_recal.csv',index=False)
    pd.set_option('display.width',200)
    print(df.pivot_table(index=['order','B','d'],columns='N',values='recall').round(2).to_string())
