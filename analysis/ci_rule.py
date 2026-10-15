"""Confidence intervals for the planning rule N95 = a + b K + c M + e/d^2 (pairs, 3 replicates).
(1) scenario bootstrap; (2) + Monte Carlo noise of each simulated recall (binomial, 30 outputs x M interactions)."""
import numpy as np, pandas as pd
P = pd.read_csv('fig/simK_pairs.csv'); rng = np.random.default_rng(2026)
def n95(N, r, t=.95):
    o = np.argsort(N); N = N[o].astype(float); r = np.maximum.accumulate(r[o])
    if r[-1] < t: return np.nan
    i = int(np.argmax(r >= t))
    if i == 0: return N[0]
    f = (t - r[i - 1]) / max(r[i] - r[i - 1], 1e-9); return float(np.exp(np.log(N[i - 1]) + f * np.log(N[i] / N[i - 1])))
G = {k: (g.N.values, g.recall.values) for k, g in P.groupby(['K', 'M', 'd'])}
keys = list(G)
def fit(tab):
    t = np.array([(K, M, d, v) for (K, M, d), v in tab if np.isfinite(v)])
    X = np.column_stack([np.ones(len(t)), t[:, 0], t[:, 1], 1 / t[:, 2] ** 2]); return np.linalg.lstsq(X, t[:, 3], rcond=None)[0]
base = fit([(k, n95(*G[k])) for k in keys]); print('point estimate', base.round(3))
out = {'scenario bootstrap': [], 'scenario + Monte Carlo': []}
for b in range(4000):
    pick = [keys[i] for i in rng.integers(0, len(keys), len(keys))]
    out['scenario bootstrap'].append(fit([(k, n95(*G[k])) for k in pick]))
    noisy = []
    for k in pick:
        N, r = G[k]; n = 30 * k[1]
        noisy.append((k, n95(N, rng.binomial(n, np.clip(r, 0, 1)) / n)))
    out['scenario + Monte Carlo'].append(fit(noisy))
rows = []
for name, v in out.items():
    v = np.array(v); lo, hi = np.percentile(v, [2.5, 97.5], axis=0)
    for j, term in enumerate(['intercept', 'K', 'M', '1/d2']): rows.append(dict(method=name, term=term, estimate=base[j], lo=lo[j], hi=hi[j]))
R = pd.DataFrame(rows); R.to_csv('fig/rule_ci.csv', index=False); print(R.round(2).to_string(index=False))
# slope in K within each (M, d) cell
S = []
for (M, d), g in P.groupby(['M', 'd']):
    t = [(K, n95(*G[(K, M, d)])) for K in sorted(g.K.unique())]; t = np.array([x for x in t if np.isfinite(x[1])])
    if len(t) >= 4: S.append(dict(M=M, d=d, slope=np.polyfit(t[:, 0], t[:, 1], 1)[0]))
S = pd.DataFrame(S); print(S.round(2).to_string(index=False)); print('median slope %.2f, range %.2f-%.2f' % (S.slope.median(), S.slope.min(), S.slope.max()))
