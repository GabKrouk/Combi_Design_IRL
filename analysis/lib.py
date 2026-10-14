"""Benchmark library: full-factorial truth, sub-designs, recovery analysis (forward selection) and metrics.
Coding: factors in +/-1 (orthogonal in the full factorial), as in a classical ANOVA."""
import itertools, numpy as np
from scipy import stats

# ------------------------------------------------------------------ model terms
def terms(K, maxorder):
    return [t for o in range(1, maxorder + 1) for t in itertools.combinations(range(K), o)]

def model_cols(A, ts):
    """A: n x K binary (0/1). Returns n x len(ts) matrix of +/-1 products."""
    S = 2 * A - 1
    return np.column_stack([np.prod(S[:, list(t)], axis=1) for t in ts]).astype(float)

def all_cells(K):
    return np.array(list(itertools.product([0, 1], repeat=K)), int)

# ------------------------------------------------------------------ dataset container
class Dataset:
    """cells: C x K binary; cell_of_obs: n_obs index into cells; Y: n_obs x G responses."""
    def __init__(self, name, factors, cells, cell_of_obs, Y, resp_names):
        self.name, self.factors, self.cells = name, factors, cells
        self.cell_of_obs, self.Y, self.resp_names = np.asarray(cell_of_obs), np.asarray(Y, float), list(resp_names)
        self.K = cells.shape[1]
        C = len(cells)
        self.cell_means = np.vstack([self.Y[self.cell_of_obs == c].mean(0) for c in range(C)])
        self.n_per_cell = np.bincount(self.cell_of_obs, minlength=C)

# ------------------------------------------------------------------ truth: saturated model
def bh(p):
    p = np.asarray(p); o = np.argsort(p, axis=None); m = p.size
    q = np.empty(m); ranked = p.ravel()[o] * m / np.arange(1, m + 1)
    q[o] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.minimum(q, 1).reshape(p.shape)

def truth_fit(ds, maxorder=None):
    """Saturated (or up to maxorder) model on all observations. Returns dict with coef, se, p, q, sigma, terms."""
    K = ds.K; ts = terms(K, maxorder or K)
    X = np.column_stack([np.ones(len(ds.cell_of_obs)), model_cols(ds.cells[ds.cell_of_obs], ts)])
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ ds.Y
    R = ds.Y - X @ B
    df = X.shape[0] - np.linalg.matrix_rank(X)
    s2 = (R ** 2).sum(0) / df
    se = np.sqrt(np.outer(np.diag(XtXi)[1:], s2))
    tval = B[1:] / se
    p = 2 * stats.t.sf(np.abs(tval), df)
    return dict(terms=ts, coef=B[1:], se=se, p=p.T, q=bh(p.T), sigma=np.sqrt(s2), df=df, intercept=B[0])

# ------------------------------------------------------------------ recovery analysis on a sub-design
def recover(ds, cell_idx, cand_order=3, alpha_main=0.05, alpha_int=0.01, max_steps=8, genes=None):
    """Main effects always in the model; interactions (order 2..cand_order) added forward while the best
    one passes F test with Bonferroni over the candidate count. Returns per gene: selected interaction
    term indices (into terms(K,cand_order)), main-effect p-values, fitted coefficient dict, residual sigma."""
    K = ds.K
    obs = np.isin(ds.cell_of_obs, cell_idx)
    A = ds.cells[ds.cell_of_obs[obs]]
    Y = ds.Y[obs] if genes is None else ds.Y[obs][:, genes]
    n, G = Y.shape
    allt = terms(K, cand_order)
    mains = [i for i, t in enumerate(allt) if len(t) == 1]
    cands = [i for i, t in enumerate(allt) if len(t) >= 2]
    Xall = model_cols(A, allt)
    B = np.column_stack([np.ones(n), Xall[:, mains]])
    Q, _ = np.linalg.qr(B)
    k0 = np.linalg.matrix_rank(B)
    Z = Xall[:, cands]
    Zp = Z - Q @ (Q.T @ Z)
    R = Y - Q @ (Q.T @ Y)
    P = len(cands)
    nrm = (Zp ** 2).sum(0)
    ok = nrm > 1e-8 * n
    sel = [[] for _ in range(G)]
    df1 = n - k0 - 1
    if df1 >= 2 and ok.any():
        S = np.where(ok[:, None], (Zp.T @ R) ** 2 / np.where(ok, nrm, 1)[:, None], -1)   # P x G
        best = S.argmax(0); bs = S.max(0); rss = (R ** 2).sum(0)
        F = bs / np.maximum(rss - bs, 1e-12) * df1
        go = stats.f.sf(F, 1, df1) * P < alpha_int
        for g in np.where(go)[0]:
            r = R[:, g].copy(); Zg = Zp.copy(); alive = ok.copy(); k = k0
            while True:
                nr = (Zg ** 2).sum(0); a = alive & (nr > 1e-8 * n)
                if not a.any() or n - k - 2 < 2 or len(sel[g]) >= max_steps: break
                sc = np.where(a, (Zg.T @ r) ** 2 / np.where(a, nr, 1), -1)
                j = int(sc.argmax()); df = n - k - 1
                Fv = sc[j] / max(r @ r - sc[j], 1e-12) * df
                if stats.f.sf(Fv, 1, df) * P >= alpha_int: break
                q = Zg[:, j] / np.sqrt(nr[j]); r = r - q * (q @ r); Zg = Zg - np.outer(q, q @ Zg)
                alive[j] = False; sel[g].append(cands[j]); k += 1
    # final fits (group genes by selected set to vectorise)
    main_p = np.ones((G, K)); coefs = [None] * G; sig = np.zeros(G)
    groups = {}
    for g in range(G): groups.setdefault(tuple(sorted(sel[g])), []).append(g)
    for sset, gs in groups.items():
        cols = mains + list(sset)
        X = np.column_stack([np.ones(n), Xall[:, cols]])
        XtXi = np.linalg.pinv(X.T @ X); Bc = XtXi @ X.T @ Y[:, gs]
        Rr = Y[:, gs] - X @ Bc; df = n - np.linalg.matrix_rank(X)
        if df < 1: s2 = np.full(len(gs), np.nan)
        else: s2 = (Rr ** 2).sum(0) / df
        se = np.sqrt(np.outer(np.diag(XtXi)[1:1 + K], s2))
        pm = 2 * stats.t.sf(np.abs(Bc[1:1 + K] / se), max(df, 1))
        for ii, g in enumerate(gs):
            main_p[g] = pm[:, ii]; sig[g] = np.sqrt(s2[ii])
            coefs[g] = (Bc[0, ii], dict(zip(cols, Bc[1:, ii])))
    return dict(sel=sel, main_p=main_p, coefs=coefs, sigma=sig, terms=allt, n=n)

def predict_cells(res, cells, K):
    """Predicted cell means for the given cells from recover() output."""
    allt = res['terms']; Xc = model_cols(cells, allt)
    out = np.zeros((len(cells), len(res['coefs'])))
    for g, (b0, d) in enumerate(res['coefs']):
        v = np.full(len(cells), b0)
        for j, b in d.items(): v += b * Xc[:, j]
        out[:, g] = v
    return out

# ------------------------------------------------------------------ truth sets and scoring
def truth_sets(tr, K, cand_order=3, q=0.01):
    """Per response: set of true mains (indices 0..K-1) and true interactions (indices in terms(K,cand_order)),
    plus count of true terms of order > cand_order."""
    allt = terms(K, cand_order); idx = {t: i for i, t in enumerate(allt)}
    sig = tr['q'] < q
    mains, ints, high = [], [], []
    for g in range(sig.shape[0]):
        m, it, h = set(), set(), 0
        for j, t in enumerate(tr['terms']):
            if not sig[g, j]: continue
            if len(t) == 1: m.add(t[0])
            elif len(t) <= cand_order: it.add(idx[t])
            else: h += 1
        mains.append(m); ints.append(it); high.append(h)
    return mains, ints, np.array(high)

def score(res, tmains, tints, alpha_main=0.05):
    G = len(tmains); out = np.zeros((G, 6))
    for g in range(G):
        foundm = set(np.where(res['main_p'][g] < alpha_main)[0])
        s = set(res['sel'][g])
        tm, ti = tmains[g], tints[g]
        main_rec = len(foundm & tm) / len(tm) if tm else np.nan
        int_rec = len(s & ti) / len(ti) if ti else np.nan
        false_int = len(s - ti)
        false_main = len(foundm - tm)
        ok = tm <= foundm and s == ti
        out[g] = [ok, main_rec, int_rec, false_int, false_main, len(s)]
    return out   # columns: success, main recall, int recall, n false ints, n false mains, n selected
