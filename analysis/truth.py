import numpy as np, pandas as pd, pickle, gzip
from lib import truth_fit, truth_sets, terms
import datasets as dsm
out={}
def stats_for(ds, tr, label, q=0.01):
    K=ds.K; tm,ti,hi=truth_sets(tr,K,cand_order=3,q=q)
    allt=terms(K,3)
    S=np.array([len(m) for m in tm]); P=np.array([sum(len(allt[j])==2 for j in s) for s in ti]); T3=np.array([sum(len(allt[j])==3 for j in s) for s in ti])
    anyeff=(S+P+T3+hi)>0
    # effect sizes relative to noise (coef in +/-1 coding = half the difference present-absent)
    sig=tr['q']<q; ratio=np.abs(tr['coef'].T)/tr['sigma'][:,None]
    order=np.array([len(t) for t in tr['terms']])
    eff={o: np.median(2*ratio[sig & (order==o)[None,:]]) if (sig & (order==o)[None,:]).any() else np.nan for o in range(1,K+1)}
    # heredity: true pair whose two factors both lack a true main effect
    pairs_tot=0; pure=0; weak=0
    for g in range(len(tm)):
        for j in ti[g]:
            t=allt[j]
            if len(t)!=2: continue
            pairs_tot+=1; n_main=sum(f in tm[g] for f in t)
            pure+= n_main==0; weak+= n_main==1
    row=dict(dataset=label,q=q,responses=len(S),with_effect=int(anyeff.sum()),
             frac_with_interaction=float(((P+T3+hi)>0)[anyeff].mean()) if anyeff.any() else np.nan,
             mean_S=float(S[anyeff].mean()), mean_pairs=float(P[anyeff].mean()), mean_triplets=float(T3[anyeff].mean()),
             mean_higher=float(hi[anyeff].mean()), max_terms=int((S+P+T3+hi).max()),
             pairs_no_main=pure/max(pairs_tot,1), pairs_one_main=weak/max(pairs_tot,1),
             eff_main=eff[1], eff_pair=eff[2], eff_trip=eff[3], sigma_median=float(np.median(tr['sigma'])))
    return row, (S,P,T3,hi)
rows=[]; dists={}
for name,f in [('Ristova',dsm.ristova),('Kellermeier_traits',dsm.kellermeier),('Kellermeier_tx',dsm.kellermeier_tx),('Tekin',dsm.tekin)]:
    ds=f(); tr=truth_fit(ds); out[name]=tr
    for q in (0.01,0.05):
        r,d=stats_for(ds,tr,name,q); rows.append(r)
        if q==0.01: dists[name]=d
    print(name,'done', flush=True)
pickle.dump((out,dists),open('prep/truth.pkl','wb'))
df=pd.DataFrame(rows); pd.set_option('display.width',250); print(df.round(3).to_string(index=False))
df.to_csv('prep/truth_summary.csv',index=False)
# sanity vs paper: Ristova probes with any term raw p<0.01 among ALL 22810 probes
with gzip.open('../data/ristova/GSE71737_series_matrix.txt.gz','rt') as h:
    for line in h:
        if line.startswith('!series_matrix_table_begin'): break
    full=pd.read_csv(h,sep='\t',index_col=0,comment='!')
from lib import Dataset, all_cells
D=pd.read_pickle('prep/ristova_design.pkl'); ci={tuple(c):i for i,c in enumerate(all_cells(5))}
dsa=Dataset('all',list('abcde'),all_cells(5),[ci[tuple(r)] for r in D.iloc[:,:5].values],np.log2(full.clip(lower=1)).T.values,full.index)
tra=truth_fit(dsa); print("Ristova all probes with any term p<0.01:",int((tra['p']<0.01).any(1).sum()),"(paper: 9507)")
