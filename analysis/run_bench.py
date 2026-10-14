"""Benchmark: sub-designs of real full factorials vs the full-factorial truth.
python run_bench.py <dataset> <out.pkl>   dataset in {ristova, tekin, kellermeier}"""
import sys, pickle, time, warnings, numpy as np
warnings.filterwarnings('ignore')
import datasets as dsm
from lib import recover, truth_fit, truth_sets, score, predict_cells, model_cols, terms, Dataset
from designs import sequences

CA3 = [3, 4, 9, 10, 15, 16, 21, 22, 27, 28]           # optimal strength-3 array for 5 factors (10 runs)
NGRID = [6, 8, 10, 12, 14, 16, 20, 24, 28, 32]
WAVES = [10, 16, 24, 32]

def simulate_like(ds, tr, tm, ti, rng):
    """Parametric copy of the data: true (FDR-significant) coefficients + Gaussian noise of the measured size."""
    allt = tr['terms']; X = model_cols(ds.cells[ds.cell_of_obs], allt)
    keep = tr['q'] < 0.01                                    # G x T
    B = np.where(keep.T, tr['coef'], 0.0)                    # T x G
    mu = tr['intercept'][None, :] + X @ B
    return Dataset(ds.name + ' (sim)', ds.factors, ds.cells, ds.cell_of_obs,
                   mu + rng.normal(size=mu.shape) * tr['sigma'][None, :], ds.resp_names)

def evaluate(ds_fit, cells_used, tm, ti, eff, truth_means, sig_cell):
    res = recover(ds_fit, cells_used)
    sc = score(res, tm, ti)
    out = dict(success=np.nanmean(sc[eff, 0]), main_rec=np.nanmean(sc[eff, 1]),
               int_rec=np.nanmean(sc[eff, 2]) if np.isfinite(sc[eff, 2]).any() else np.nan,
               false_int_eff=np.mean(sc[eff, 3] > 0), false_int_null=np.mean(sc[~eff, 3] > 0) if (~eff).any() else np.nan,
               n_sel=np.mean(sc[:, 5]))
    # success split by complexity of the truth
    nint = np.array([len(x) for x in ti])
    out['success_noint'] = np.nanmean(sc[eff & (nint == 0), 0]) if (eff & (nint == 0)).any() else np.nan
    out['success_int'] = np.nanmean(sc[eff & (nint > 0), 0]) if (eff & (nint > 0)).any() else np.nan
    held = np.setdiff1d(np.arange(len(ds_fit.cells)), cells_used)
    if len(held):
        pred = predict_cells(res, ds_fit.cells[held], ds_fit.K)
        err = (pred - truth_means[held]) / sig_cell[None, :]
        out['pred_ratio'] = float(np.median(np.sqrt((err[:, eff] ** 2).mean(0))))
    else: out['pred_ratio'] = np.nan
    return out, res, sc

def stopping(ds_fit, seq, tm, ti, eff, thr=1.5):
    """Wave rule per response: fit on waves so far, predict next wave; stop (after measuring it) when the
    prediction error is within thr x expected noise of a cell mean. Returns success and mixtures used."""
    G = ds_fit.Y.shape[1]; stopped = np.zeros(G, bool); used = np.full(G, WAVES[-1]); succ = np.zeros(G)
    reps = np.median(ds_fit.n_per_cell)
    for w in range(len(WAVES) - 1):
        cur, nxt = seq[:WAVES[w]], seq[WAVES[w]:WAVES[w + 1]]
        res = recover(ds_fit, cur)
        pred = predict_cells(res, ds_fit.cells[nxt], ds_fit.K)
        obs = ds_fit.cell_means[nxt]
        ratio = np.sqrt(((pred - obs) ** 2).mean(0)) / (res['sigma'] * np.sqrt(1 / reps + 1 / reps))
        stop_now = (~stopped) & (ratio < thr)
        if stop_now.any():
            res2 = recover(ds_fit, seq[:WAVES[w + 1]]); sc = score(res2, tm, ti)
            succ[stop_now] = sc[stop_now, 0]; used[stop_now] = WAVES[w + 1]; stopped |= stop_now
    if (~stopped).any():
        sc = score(recover(ds_fit, seq), tm, ti); succ[~stopped] = sc[~stopped, 0]
    return dict(stop_success=float(succ[eff].mean()), stop_N=float(used[eff].mean()),
                stop_trap=float(np.mean((used[eff] < WAVES[-1]) & (succ[eff] == 0))))

def lasso_eval(ds_fit, cells_used, tm, ti, idx):
    from sklearn.linear_model import LassoLarsIC
    K = ds_fit.K; allt = terms(K, 3); mains = [i for i, t in enumerate(allt) if len(t) == 1]
    cands = [i for i, t in enumerate(allt) if len(t) >= 2]
    obs = np.isin(ds_fit.cell_of_obs, cells_used)
    X = model_cols(ds_fit.cells[ds_fit.cell_of_obs[obs]], allt); n = X.shape[0]
    B = np.column_stack([np.ones(n), X[:, mains]]); Q, _ = np.linalg.qr(B)
    Z = X[:, cands] - Q @ (Q.T @ X[:, cands]); nz = np.linalg.norm(Z, axis=0); okc = nz > 1e-8
    succ = []
    for g in idx:
        y = ds_fit.Y[obs, g]; r = y - Q @ (Q.T @ y)
        try:
            m = LassoLarsIC(criterion='bic', fit_intercept=False).fit(Z[:, okc] / nz[okc], r)
            sel = {cands[j] for j, c in zip(np.where(okc)[0], m.coef_) if abs(c) > 1e-10}
        except Exception: sel = set()
        succ.append(sel == ti[g])
    return float(np.mean(succ))

def full_design_sets(ds_truth):
    """Conclusions of the same analysis applied to the complete factorial (second reference)."""
    res = recover(ds_truth, np.arange(len(ds_truth.cells)))
    tm = [set(np.where(res['main_p'][g] < 0.05)[0]) for g in range(len(res['sel']))]
    ti = [set(x) for x in res['sel']]
    return tm, ti

def run_one(ds_fit, tr_truth, truth_means, rng, draws, do_sim=True, do_lasso=False, ds_truth=None):
    tm, ti, hi = truth_sets(tr_truth, ds_fit.K)
    tmF, tiF = full_design_sets(ds_truth if ds_truth is not None else ds_fit)
    effF = np.array([len(a) + len(b) > 0 for a, b in zip(tmF, tiF)])
    eff = np.array([len(a) + len(b) + h > 0 for a, b, h in zip(tm, ti, hi)])
    sig_cell = tr_truth['sigma'] / np.sqrt(np.median(ds_fit.n_per_cell))
    rows = []
    lidx = None
    if do_lasso:
        e = np.where(eff)[0]; lidx = rng.choice(e, min(400, len(e)), replace=False)
        tm_l = [tm[i] for i in range(len(tm))]
    for d in range(draws):
        seqs = sequences(ds_fit.K, rng, CA3)
        sim = simulate_like(ds_fit, tr_truth, tm, ti, rng) if do_sim else None
        for dname, seq in seqs.items():
            for N in NGRID:
                cells = np.array(seq[:N])
                o, _, _ = evaluate(ds_fit, cells, tm, ti, eff, truth_means, sig_cell)
                o.update(design=dname, N=N, draw=d, kind='real')
                oF, _, _ = evaluate(ds_fit, cells, tmF, tiF, effF, truth_means, sig_cell)
                o['agree_full'] = oF['success']; o['agree_full_int'] = oF['success_int']
                if do_lasso and d < 5 and dname != 'one_at_a_time':
                    o['lasso_int_exact'] = lasso_eval(ds_fit, cells, tm, ti, lidx)
                    rs = recover(ds_fit, cells, genes=lidx)
                    o['fs_int_exact'] = float(np.mean([set(rs['sel'][k]) == ti[g] for k, g in enumerate(lidx)]))
                rows.append(o)
                if do_sim:
                    o2, _, _ = evaluate(sim, cells, tm, ti, eff, truth_means, sig_cell)
                    o2.update(design=dname, N=N, draw=d, kind='simulated'); rows.append(o2)
            if dname == 'combinatorial':
                s = stopping(ds_fit, seq, tm, ti, eff); s.update(design=dname, N=-1, draw=d, kind='stop'); rows.append(s)
    return rows

if __name__ == '__main__':
    name, outp = sys.argv[1], sys.argv[2]
    rng = np.random.default_rng(12345)
    t0 = time.time(); rows = []
    if name in ('ristova', 'tekin'):
        ds = dsm.ristova() if name == 'ristova' else dsm.tekin()
        if name == 'tekin':          # 2,000 blocks drawn at random to keep the run time reasonable
            sub = np.sort(np.random.default_rng(7).choice(ds.Y.shape[1], 2000, replace=False))
            ds = Dataset(ds.name, ds.factors, ds.cells, ds.cell_of_obs, ds.Y[:, sub], [ds.resp_names[i] for i in sub])
        tr = truth_fit(ds)
        rows = run_one(ds, tr, ds.cell_means, rng, draws=30, do_lasso=True)
    elif name == 'kellermeier':
        for split in range(20):          # truth from half A, designs analysed on half B (independent plants)
            A = dsm.kellermeier(seed=split, half=0); Bh = dsm.kellermeier(seed=split, half=1)
            trA = truth_fit(A)
            r = run_one(Bh, trA, A.cell_means, rng, draws=5, do_sim=True, ds_truth=A)
            for x in r: x['split'] = split
            rows += r
            print('split', split, '%.0fs' % (time.time() - t0), flush=True)
    pickle.dump(rows, open(outp, 'wb'))
    print(name, 'rows', len(rows), '%.0fs' % (time.time() - t0))
