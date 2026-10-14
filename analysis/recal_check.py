"""Two-level recalibrated simulation (true effects drawn from the fitted mixtures, 'truth' re-estimated from the
simulated full factorial exactly as in the real benchmark) and comparison with the real recall by order."""
import numpy as np, pandas as pd, pickle, warnings; warnings.filterwarnings('ignore')
import datasets as dsm
from lib import Dataset, truth_fit, truth_sets, recover, terms, model_cols
from designs import sequences
from run_bench import CA3
fits=pickle.load(open('prep/eb_fits.pkl','rb'))
allt5=terms(5,5); order5=np.array([len(t) for t in allt5]); allt3=terms(5,3); order3=np.array([len(t) for t in allt3])
def draw_theta(name, G, rng):
    th=np.zeros((len(allt5),G))
    for o in range(1,6):
        lab={1:'simples',2:'paires',3:'triplets'}.get(o,'ordre 4-5')
        if (name,lab) not in fits: continue
        p,t=fits[(name,lab)]; idx=np.where(order5==o)[0]
        comp=rng.choice(3,size=(len(idx),G),p=p/p.sum())
        th[idx]=np.where(comp==0,0,rng.normal(size=(len(idx),G))*t[comp])
    return th
def sim_ds(template, th, rng, noise=None):
    X=model_cols(template.cells[template.cell_of_obs],allt5)
    E=rng.normal(size=(X.shape[0],th.shape[1])) if noise is None else rng.choice(noise,size=(X.shape[0],th.shape[1]))
    Y=X@th+E
    return Dataset('sim',template.factors,template.cells,template.cell_of_obs,Y,range(th.shape[1]))
def recall(res,ti,idx):
    out={}
    for o in (2,3):
        hit=tot=0
        for g in idx:
            tg={j for j in ti[g] if order3[j]==o}; sg={j for j in res['sel'][g] if order3[j]==o}; hit+=len(tg&sg); tot+=len(tg)
        out[o]=hit/tot if tot else np.nan
    return out
rows=[]
NOISE={}
def run(name, tmpl, G, draws, rng, halves=False, tmplB=None):
    th=draw_theta(name,G,rng)
    nz=NOISE.get(name)
    A=sim_ds(tmpl,th,rng,nz); B=sim_ds(tmplB,th,rng,nz) if halves else A
    tr=truth_fit(A); tm,ti,hi=truth_sets(tr,5); idx=[g for g in range(G) if ti[g]]
    for d in range(draws):
        seqs=sequences(5,rng,CA3)
        for dname in ('combinatorial','random','one_at_a_time'):
            for N in (10,16,24,32):
                r=recall(recover(B,np.array(seqs[dname][:N])),ti,idx)
                for o in (2,3): rows.append(dict(jeu=name,plan=dname,N=N,ordre=o,rappel=r[o]))
import sys
rng=np.random.default_rng(11)
if len(sys.argv)>1 and sys.argv[1]=='empirical':
    def resid(ds):
        tr=truth_fit(ds); X=np.column_stack([np.ones(len(ds.cell_of_obs)),model_cols(ds.cells[ds.cell_of_obs],tr['terms'])])
        R=ds.Y-X@np.vstack([tr['intercept'],tr['coef']]); n,p=X.shape; ok=tr['sigma']>0
        z=(R[:,ok]/tr['sigma'][ok])*np.sqrt(n/(n-p)); return z[np.abs(z)>1e-9].ravel()
    NOISE['Ristova']=resid(dsm.ristova()); NOISE['Tekin']=resid(dsm.tekin()); NOISE['Kellermeier traits']=resid(dsm.kellermeier())
    for k,v in NOISE.items(): print(k,'kurtosis excess %.2f'%(((v-v.mean())**4).mean()/v.var()**2-3))
run('Ristova',dsm.ristova(),3000,6,rng)
t=dsm.tekin(); run('Tekin',t,2000,6,rng)
for sp in range(5):
    run('Kellermeier traits',dsm.kellermeier(seed=sp,half=0),300,2,rng,halves=True,tmplB=dsm.kellermeier(seed=sp,half=1))
S=pd.DataFrame(rows).groupby(['jeu','ordre','plan','N'])['rappel'].mean().rename('recalé')
R=pd.read_csv('fig/recall_by_order.csv'); R['jeu']=R['jeu'].replace({'Kellermeier':'Kellermeier traits'})
real=R[R.type=='réel'].groupby(['jeu','ordre','plan','N'])['rappel'].mean().rename('réel')
old=R[R.type=='simulé'].groupby(['jeu','ordre','plan','N'])['rappel'].mean().rename('ancien_sim')
T=pd.concat([real,old,S],axis=1).reset_index(); T.to_csv('fig/recalibration_check%s.csv'%('_emp' if len(sys.argv)>1 else ''),index=False)
pd.set_option('display.width',200); print(T[T.plan=='combinatorial'].round(2).to_string(index=False))
T=T[~((T.jeu=='Ristova')&(T.ordre==3))]
for c in ('ancien_sim','recalé'):
    e=(T[c]-T['réel']).dropna(); print(c,'erreur moyenne %.3f, |erreur| moyenne %.3f'%(e.mean(),e.abs().mean()))
