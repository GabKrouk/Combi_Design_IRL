"""Build the benchmark datasets (Dataset objects) from the prepared files."""
import pickle, re, numpy as np, pandas as pd
from lib import Dataset, all_cells

def cell_index(cells):
    return {tuple(c): i for i, c in enumerate(cells)}

def ristova():
    X = pd.read_pickle('prep/ristova_log2.pkl'); D = pd.read_pickle('prep/ristova_design.pkl')
    fac = ['NO3', 'NH4', 'IAA', 'CK', 'ABA']; cells = all_cells(5); ci = cell_index(cells)
    cob = [ci[tuple(r)] for r in D[fac].values]
    return Dataset('Ristova 2016 (transcriptome)', fac, cells, cob, X.T.values, X.index)

def kellermeier(seed=0, half=None):
    d = pd.read_pickle('prep/kellermeier_ds1.pkl'); d = d[d['Age'] == 10].copy()
    fac = ['N', 'P', 'K', 'S', 'L']; cells = all_cells(5); ci = cell_index(cells)
    traits = ['Total Root Size', 'MR length', 'Basal', 'Branched', 'Apical', 'MR angle', 'Lateral Root Size',
              '1st order LR #', '2nd order LR #', 'LRP 0.25', 'LRP 0.50', 'LRP 0.75', 'LRP 1.00',
              'LR density / MR', 'LR density / BZ']
    Y = d[traits].astype(float).values
    Y = np.where(np.array([t != 'MR angle' for t in traits]), np.log1p(np.clip(Y, 0, None)), Y)
    cob = np.array([ci[tuple(r)] for r in d[fac].values])
    ok = ~np.isnan(Y).any(1); Y, cob = Y[ok], cob[ok]
    if half is not None:          # split plants within each condition into two halves
        rng = np.random.default_rng(seed); keep = np.zeros(len(cob), bool)
        for c in range(32):
            ii = np.where(cob == c)[0]; rng.shuffle(ii); h = ii[: len(ii) // 2] if half == 0 else ii[len(ii) // 2:]
            keep[h] = True
        Y, cob = Y[keep], cob[keep]
    return Dataset('Kellermeier 2014 (root traits)', fac, cells, cob, Y, traits)

def kellermeier_tx():
    L = pd.read_pickle('prep/kellermeier_ds4_log2.pkl'); D = pd.read_pickle('prep/kellermeier_ds4_design.pkl')
    fac = ['N', 'P', 'K']; cells = all_cells(3); ci = cell_index(cells)
    cob = [ci[tuple(r)] for r in D[fac].values]
    return Dataset('Kellermeier 2014 (transcriptome NPK)', fac, cells, cob, L.values.T, L.index)

def tekin():
    blocks = pickle.load(open('prep/tekin_blocks.pkl', 'rb'))
    cells = all_cells(5); ci = cell_index(cells)
    names0 = blocks[0][1]
    Ys, codes = [], []
    for code, names, vals in blocks:
        drugs = re.findall(r'[A-Z]{3}', code)
        cols = []
        for nm in names:
            v = [0] * 5
            for dname in nm.split('+'): v[drugs.index(dname)] = 1
            cols.append(ci[tuple(v)])
        y = np.log2(np.clip(vals, 1, None) / 100.0)          # 3 x 31, relative growth on log scale
        full = np.zeros((3, 32)); full[:, cols] = y            # control cell (index 0) stays log2(1) = 0
        Ys.append(full.ravel()); codes.append(code)
    cob = np.tile(np.arange(32), 3)
    return Dataset('Tekin 2018 (E. coli, 5 antibiotics)', ['d1', 'd2', 'd3', 'd4', 'd5'], cells, cob, np.array(Ys).T, codes)

if __name__ == '__main__':
    for f in (ristova, lambda: kellermeier(), kellermeier_tx, tekin):
        ds = f(); print(ds.name, ds.Y.shape, 'obs per cell', np.unique(ds.n_per_cell))
