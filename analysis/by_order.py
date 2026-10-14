import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
import datasets as dsm
from lib import truth_fit, truth_sets, recover, terms, Dataset
from designs import sequences
from run_bench import CA3, simulate_like
allt=terms(5,3); order=np.array([len(t) for t in allt])
def recall_by_order(res, ti, eff_idx):
    out={}
    for o in (2,3):
        hit=tot=0; false=0
        for g in eff_idx:
            tg={j for j in ti[g] if order[j]==o}; sg={j for j in res['sel'][g] if order[j]==o}
            hit+=len(tg&sg); tot+=len(tg); false+=len(sg-tg)
        out[o]=(hit/tot if tot else np.nan, tot, false)
    return out
rows=[]
def run(name, ds_fit, tr, draws, ds_truth=None, rng=np.random.default_rng(3)):
    tm,ti,hi=truth_sets(tr,5); idx=[g for g in range(len(ti)) if ti[g]]
    for d in range(draws):
        seqs=sequences(5,rng,CA3); sim=simulate_like(ds_fit,tr,tm,ti,rng)
        for dname in ('combinatorial','random','one_at_a_time'):
            for N in (10,16,24,32):
                cells=np.array(seqs[dname][:N])
                for kind,D in (('réel',ds_fit),('simulé',sim)):
                    r=recall_by_order(recover(D,cells),ti,idx)
                    for o in (2,3): rows.append(dict(jeu=name,plan=dname,N=N,type=kind,ordre=o,rappel=r[o][0],n_vrais=r[o][1],faux=r[o][2]))
ds=dsm.ristova(); run('Ristova',ds,truth_fit(ds),10)
ds=dsm.tekin(); sub=np.sort(np.random.default_rng(7).choice(ds.Y.shape[1],2000,replace=False))
ds=Dataset(ds.name,ds.factors,ds.cells,ds.cell_of_obs,ds.Y[:,sub],[ds.resp_names[i] for i in sub]); run('Tekin',ds,truth_fit(ds),10)
for sp in range(10):
    A=dsm.kellermeier(seed=sp,half=0); B=dsm.kellermeier(seed=sp,half=1); run('Kellermeier',B,truth_fit(A),2)
df=pd.DataFrame(rows); df.to_csv('fig/recall_by_order.csv',index=False)
pd.set_option('display.width',200)
t=df.groupby(['jeu','ordre','type','plan','N'])['rappel'].mean().unstack('N').round(2)
print(t.to_string())
print(df[(df.type=='réel')&(df.plan=='combinatorial')&(df.N==32)].groupby(['jeu','ordre'])[['n_vrais']].mean())
