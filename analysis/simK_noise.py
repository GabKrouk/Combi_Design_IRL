"""K sweep (pairs, M=5, d=1) with the empirical, heavy-tailed noise of each real dataset (standardized residuals of
the complete factorial fits) instead of Gaussian noise. Gives a dataset-based interval for the K coefficient."""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
import datasets as dsm
from lib import Dataset, recover, terms, model_cols, truth_fit

def residuals(ds):
    tr = truth_fit(ds); X = np.column_stack([np.ones(len(ds.cell_of_obs)), model_cols(ds.cells[ds.cell_of_obs], tr['terms'])])
    R = ds.Y - X @ np.vstack([tr['intercept'], tr['coef']]); n, p = X.shape; ok = tr['sigma'] > 0
    z = (R[:, ok] / tr['sigma'][ok]) * np.sqrt(n / (n - p)); z = z[np.abs(z) > 1e-9].ravel()
    return z / z.std()

def one(args):
    K, M, d, B, N, seed, G, noise = args
    rng = np.random.default_rng(seed); R = 3
    A = np.zeros((N, K), int)
    for j in range(K): A[rng.permutation(N)[:N // 2], j] = 1
    allt = terms(K, 2); o = np.array([len(t) for t in allt]); X = model_cols(np.repeat(A, R, 0), allt)
    cand = np.where(o == 2)[0]; th = np.zeros((len(allt), G))
    for g in range(G):
        m = rng.choice(K, K // 2, replace=False); th[m, g] = rng.normal(0, 1.1, K // 2)
        tg = rng.choice(cand, M, replace=False); th[tg, g] = rng.choice([-1, 1], M) * d / 2
        bg = np.setdiff1d(cand, tg); th[bg, g] += rng.normal(0, np.sqrt(B / len(bg)), len(bg))
    E = rng.normal(size=(X.shape[0], G)) if noise is None else rng.choice(NOISE[noise], size=(X.shape[0], G))
    Y = X @ th + E
    ds = Dataset('s', list(range(K)), A, np.repeat(np.arange(N), R), Y, range(G))
    res = recover(ds, np.arange(N), cand_order=2, max_steps=M + 8)
    rec = [len(set(res['sel'][g]) & set(np.where((np.abs(th[:, g]) >= d / 2 - 1e-12) & (o == 2))[0])) / M for g in range(G)]
    return float(np.mean(rec))

NOISE = {}
def scen(a):
    K, noise, seed = a; rows = []; N = K + 12; hits = 0
    while N <= 15 * K + 60:
        r = one((K, 5, 1.0, 0.1, int(N), seed + int(N), 30, noise)); rows.append(dict(K=K, noise=noise or 'gaussian', N=int(N), recall=r))
        hits = hits + 1 if r >= .97 else 0
        if hits >= 2: break
        N = max(N + 2, int(round(N * 1.2)))
    return rows

if __name__ == '__main__':
    NOISE['ristova'] = residuals(dsm.ristova()); NOISE['kellermeier'] = residuals(dsm.kellermeier()); NOISE['tekin'] = residuals(dsm.tekin())
    for k, v in NOISE.items(): print(k, 'excess kurtosis %.1f' % (((v - v.mean()) ** 4).mean() / v.var() ** 2 - 3), flush=True)
    S = [(K, nz, 1000 * i + 77) for i, (K, nz) in enumerate((K, nz) for K in (100, 75, 50, 30, 18, 10) for nz in (None, 'ristova', 'kellermeier', 'tekin'))]
    rows = []
    with Pool(2) as p:
        for r in p.imap_unordered(scen, S):
            rows += r; pd.DataFrame(rows).to_csv('fig/simK_noise.csv', index=False); print(len(rows), flush=True)
    print('done')
