"""Same analysis applied to the exhaustive pairwise design (control, every molecule alone, every pair), 3 replicates."""
import numpy as np, pandas as pd, sys, time
from lib import Dataset, recover, terms, model_cols
R = 3
def one_ex(K, M, d, B, seed, G=30):
    rng = np.random.default_rng(seed)
    A = [np.zeros(K, int)] + [np.eye(K, dtype=int)[i] for i in range(K)]
    for i in range(K):
        for j in range(i + 1, K):
            v = np.zeros(K, int); v[[i, j]] = 1; A.append(v)
    A = np.array(A); N = len(A)
    allt = terms(K, 2); o = np.array([len(t) for t in allt]); X = model_cols(np.repeat(A, R, 0), allt)
    cand = np.where(o == 2)[0]; th = np.zeros((len(allt), G)); S = K // 2
    for g in range(G):
        m = rng.choice(K, S, replace=False); th[m, g] = rng.normal(0, 1.1, S)
        tgt = rng.choice(cand, M, replace=False); th[tgt, g] = rng.choice([-1, 1], M) * d / 2
        bg = np.setdiff1d(cand, tgt); th[bg, g] += rng.normal(0, np.sqrt(B / len(bg)), len(bg))
    Y = X @ th + rng.normal(size=(X.shape[0], G))
    ds = Dataset('ex', list(range(K)), A, np.repeat(np.arange(N), R), Y, range(G))
    res = recover(ds, np.arange(N), cand_order=2, max_steps=M + 8)
    rec = []; fdr = []
    for g in range(G):
        tgt = set(np.where((np.abs(th[:, g]) >= d / 2 - 1e-12) & (o == 2))[0]); s = set(res['sel'][g]); hit = len(s & tgt)
        rec.append(hit / M); fdr.append((len(s) - hit) / len(s) if s else 0)
    return dict(K=K, M=M, d=d, N=N, recall=np.mean(rec), fdr=np.mean(fdr))
if __name__ == '__main__':
    rows = []
    for K in (10, 18, 30, 50):
        for d in (0.5, 1.0, 2.0):
            t = time.time(); r = one_ex(K, 5, d, 0.1, 31 * K + int(10 * d)); rows.append(r)
            pd.DataFrame(rows).to_csv('fig/simK_exhaustive.csv', index=False); print(r, '%.0fs' % (time.time() - t), flush=True)
