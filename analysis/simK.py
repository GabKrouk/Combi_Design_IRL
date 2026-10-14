"""Planning simulations for K = 10..100 compounds (random balanced designs, 3 replicates).
Effects in +/-1 coding; d = difference present/absent of a target interaction, in replicate-SD units.
Mains: S = K/2 active compounds, coefficients ~ N(0, 1.1). Background: small interactions with total variance B."""
import sys, time, numpy as np, pandas as pd, warnings, itertools; warnings.filterwarnings('ignore')
from multiprocessing import Pool
from lib import Dataset, recover, terms, model_cols
R = 3
def one(args):
    K, order, M, d, B, N, seed, G = args
    rng = np.random.default_rng(seed)
    A = np.zeros((N, K), int)
    for j in range(K): A[rng.permutation(N)[:N // 2], j] = 1
    allt = terms(K, order); o = np.array([len(t) for t in allt])
    X = model_cols(np.repeat(A, R, 0), allt)
    cand = np.where(o == order)[0]; bgpool = np.where(o >= 2)[0]
    th = np.zeros((len(allt), G)); S = K // 2
    for g in range(G):
        m = rng.choice(K, S, replace=False); th[m, g] = rng.normal(0, 1.1, S)
        tgt = rng.choice(cand, M, replace=False); th[tgt, g] = rng.choice([-1, 1], M) * d / 2
        if B > 0:
            bg = np.setdiff1d(bgpool, tgt); th[bg, g] += rng.normal(0, np.sqrt(B / len(bg)), len(bg))
    Y = X @ th + rng.normal(size=(X.shape[0], G))
    ds = Dataset('s', list(range(K)), A, np.repeat(np.arange(N), R), Y, range(G))
    t0 = time.time()
    res = recover(ds, np.arange(N), cand_order=order, max_steps=M + 8)
    rec = []; fdr = []; mainrec = []
    for g in range(G):
        tgt = set(np.where((np.abs(th[:, g]) >= d / 2 - 1e-12) & (o == order))[0])
        s = set(res['sel'][g]); hit = len(s & tgt); rec.append(hit / M); fdr.append((len(s) - hit) / len(s) if s else 0)
        act = np.where(np.abs(th[:K, g]) > 0.25)[0]
        mainrec.append(np.mean(res['main_p'][g][act] < 0.05) if len(act) else np.nan)
    return dict(K=K, order=order, M=M, d=d, B=B, N=N, recall=np.mean(rec), fdr=np.mean(fdr), main_recall=np.nanmean(mainrec), sec=time.time() - t0)
def ngrid(K, M, d, order):
    base = K + 12 + 2 * M
    f = [1.0, 1.25, 1.6, 2.0, 2.6, 3.4, 4.5, 6.0]
    return sorted({int(round(base * x)) for x in f})
if __name__ == '__main__':
    part = sys.argv[1]; out = sys.argv[2]; jobs = []; s = 0
    if part == 'pairs':
        for K in (10, 18, 30, 50, 75, 100):
            for M in (1, 3, 5, 10):
                for d in (0.5, 1.0, 2.0):
                    for N in ngrid(K, M, d, 2): s += 1; jobs.append((K, 2, M, d, 0.1, N, s, 40))
    elif part == 'bg':
        for K in (18, 50, 100):
            for B in (0.0, 0.45):
                for d in (0.5, 1.0, 2.0):
                    for N in ngrid(K, 5, d, 2): s += 1; jobs.append((K, 2, 5, d, B, N, 10000 + s, 40))
    elif part == 'trip':
        for K in (10, 18, 30):
            for M in (1, 3, 5):
                for d in (1.0, 2.0):
                    for N in ngrid(K, M, d, 3): s += 1; jobs.append((K, 3, M, d, 0.1, N, 20000 + s, 30))
    jobs.sort(key=lambda j: -j[0] * j[5])   # big first
    rows = []
    with Pool(int(sys.argv[3]) if len(sys.argv) > 3 else 1) as p:
        for r in p.imap_unordered(one, jobs):
            rows.append(r); pd.DataFrame(rows).to_csv(out, index=False)
    print('done', len(rows))
