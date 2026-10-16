"""Replication versus diversification for K = 18 (3 replicates per mixture vs one measurement per mixture)."""
"""Wave plan for 18 molecules and in-silico evaluation at each wave cut."""
import sys, itertools, numpy as np, pandas as pd
sys.path.insert(0, '.')
from lib import Dataset, recover, terms, model_cols

K = 18; R = 3
CUTS = [24, 48, 88, 152, 240]
PAIRS = list(itertools.combinations(range(K), 2))

def imbalance(A):
    """Sum over pairs of squared deviation of the 4 presence/absence counts from N/4, plus column balance."""
    n = len(A); s = 0.0
    for i, j in PAIRS:
        c = np.bincount(A[:, i] * 2 + A[:, j], minlength=4); s += ((c - n / 4) ** 2).sum()
    return s + 5 * ((A.sum(0) - n / 2) ** 2).sum()

def augment(A, m, rng, cand=600):
    for _ in range(m):
        best = None
        for _ in range(cand):
            r = (rng.random(K) < 0.5).astype(int)
            if r.sum() < 6 or r.sum() > 12: continue
            if any((A == r).all(1)): continue
            v = imbalance(np.vstack([A, r]))
            if best is None or v < best[0]: best = (v, r)
        A = np.vstack([A, best[1]])
    return A

def build(seed=1):
    rng = np.random.default_rng(seed)
    ca = np.load('fig/ca3_K18.npy')
    A = np.vstack([ca, np.zeros(K, int)])            # wave 1: strength-3 array + control (no VOC)
    A = augment(A, 24 - len(A), rng)                  # -> 24
    A = np.vstack([A, 1 - ca])                        # wave 2: foldover of the array (+20)
    A = augment(A, 48 - len(A), rng)                  # -> 48
    for c in CUTS[2:]: A = augment(A, c - len(A), rng)
    return A

def evaluate(A, order, M, d, B, reps=8, G=40, seed=0):
    rng = np.random.default_rng(seed); out = {c: dict(rec=[], main=[], fdr=[]) for c in CUTS}
    allt = terms(K, order); o = np.array([len(t) for t in allt]); cand = np.where(o == order)[0]; bgpool = np.where(o >= 2)[0]
    for rep in range(reps):
        perm = rng.permutation(K); D = A[:, perm]
        X = model_cols(np.repeat(D, R, 0), allt)
        th = np.zeros((len(allt), G))
        for g in range(G):
            m = rng.choice(K, 14, replace=False); th[m, g] = rng.normal(0, 1.1, 14)
            tg = rng.choice(cand, M, replace=False); th[tg, g] = rng.choice([-1, 1], M) * d / 2
            bg = np.setdiff1d(bgpool, tg); th[bg, g] += rng.normal(0, np.sqrt(B / len(bg)), len(bg))
        Y = X @ th + rng.normal(size=(X.shape[0], G))
        ds = Dataset('w', list(range(K)), D, np.repeat(np.arange(len(D)), R), Y, range(G))
        for c in CUTS:
            res = recover(ds, np.arange(c), cand_order=order, max_steps=M + 8)
            for g in range(G):
                tg = set(np.where((np.abs(th[:, g]) >= d / 2 - 1e-12) & (o == order))[0]); s = set(res['sel'][g]); h = len(s & tg)
                out[c]['rec'].append(h / M); out[c]['fdr'].append((len(s) - h) / len(s) if s else 0)
                act = np.where(np.abs(th[:K, g]) > 0.25)[0]; out[c]['main'].append(np.mean(res['main_p'][g][act] < 0.05))
    return [dict(order=order, M=M, d=d, B=B, N=c, recall=np.mean(v['rec']), main=np.mean(v['main']), fdr=np.mean(v['fdr'])) for c, v in out.items()]

if __name__ == '__main__':
    A = build(); np.save('fig/plan_K18_240.npy', A)
    print('rows', len(A), 'per-molecule presence', A.sum(0).min(), A.sum(0).max(), 'molecules per mixture', A.sum(1).min(), A.sum(1).max())
    for c in CUTS:
        cnt = min(np.bincount(A[:c, i] * 2 + A[:c, j], minlength=4).min() for i, j in PAIRS); print(c, 'min pair-state count', cnt)
    rows = []
    for order, M, d, B in [(2, 5, 0.5, 0.1), (2, 5, 1.0, 0.1), (2, 5, 2.0, 0.1), (2, 5, 1.0, 0.45), (3, 3, 1.0, 0.1), (3, 3, 2.0, 0.1)]:
        r = evaluate(A, order, M, d, B, reps=6 if order == 3 else 8); rows += r
        print(order, M, d, B, [(x['N'], round(x['recall'], 2), round(x['main'], 2)) for x in r], flush=True)
    pd.DataFrame(rows).to_csv('fig/rep_sim_R3.csv', index=False)
