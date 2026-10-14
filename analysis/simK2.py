"""Adaptive N search per scenario (pairs): geometric N grid until recall >= 0.97 twice."""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
from simK import one


def scen(a):
    K, M, d, B, seed = a
    rows = []; N = K + 2 + 2 * M; hits = 0
    while N <= 15 * K + 60:
        r = one((K, 2, M, d, B, int(N), seed + int(N), 30)); rows.append(r)
        hits = hits + 1 if r['recall'] >= 0.97 else 0
        if hits >= 2: break
        N = max(N + 2, int(round(N * 1.3)))
    return rows


if __name__ == '__main__':
    part, out = sys.argv[1], sys.argv[2]
    if part == 'pairs':
        combos = [(K, M, d) for K in (10, 18, 30, 50, 75, 100) for M in (1, 3, 5, 10) for d in (0.5, 1.0, 2.0)]
        S = [(K, M, d, 0.1, 1000 * i) for i, (K, M, d) in enumerate(combos)]
    else:
        combos = [(K, B, d) for K in (18, 50, 100) for B in (0.0, 0.45) for d in (0.5, 1.0, 2.0)]
        S = [(K, 5, d, B, 7 * 10**6 + 1000 * i) for i, (K, B, d) in enumerate(combos)]
    S.sort(key=lambda s: -s[0]); rows = []
    with Pool(2) as p:
        for r in p.imap_unordered(scen, S):
            rows += r; pd.DataFrame(rows).to_csv(out, index=False); print(len(rows), flush=True)
    print('done')
