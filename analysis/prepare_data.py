"""Turn the raw public downloads (../data) into the tables used by datasets.py (prep/*.pkl).

    python prepare_data.py                 # build everything into prep/
    python prepare_data.py --check REF     # also compare with a reference prep folder
"""
import gzip, os, re, sys, pickle, subprocess, io
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, 'prep')
os.makedirs(OUT, exist_ok=True)


# ----------------------------------------------------------------- Ristova et al. 2016 (GSE71737)
def ristova():
    """MAS5 signals from the GEO series matrix; probes (Affymetrix controls excluded) with signal > 100 in at
    least half of the 64 arrays, log2. Design from sample titles 'NO-NH-IA-CK-AB a-b-c-d-e-repN' (a value of 0.5 counts as present)."""
    path = os.path.join(DATA, 'ristova', 'GSE71737_series_matrix.txt.gz')
    titles = None
    with gzip.open(path, 'rt') as h:
        for line in h:
            if line.startswith('!Sample_title'):
                titles = [t.strip('"') for t in line.rstrip('\n').split('\t')[1:]]
            if line.startswith('!series_matrix_table_begin'):
                break
        X = pd.read_csv(h, sep='\t', index_col=0, comment='!')
    X.columns = [c.strip('"') for c in X.columns]
    rows = []
    for t in titles:
        m = re.search(r'([\d.]+)-([\d.]+)-([\d.]+)-([\d.]+)-([\d.]+)-rep(\d+)', t)
        lev = [int(float(v) > 0) for v in m.groups()[:5]]
        rows.append(lev + [int(m.group(6))])
    D = pd.DataFrame(rows, index=X.columns, columns=['NO3', 'NH4', 'IAA', 'CK', 'ABA', 'rep'])
    keep = ((X > 100).mean(axis=1) >= 0.5) & ~X.index.str.startswith('AFFX')   # drop Affymetrix controls
    L = np.log2(X[keep].clip(lower=1))
    L.to_pickle(os.path.join(OUT, 'ristova_log2.pkl')); D.to_pickle(os.path.join(OUT, 'ristova_design.pkl'))
    print('Ristova:', L.shape[0], 'probes x', L.shape[1], 'arrays')


# ----------------------------------------------------------------- Kellermeier et al. 2014 (Plant Cell)
def kellermeier():
    """DS1: raw root-architecture measurements (header on the 3rd row, units on the 4th).
    Factor coding: upper-case letter in the medium = high supply (N, P, K, S); long day = light (L = 1).
    DS4: absolute expression (2^3 N/P/K x 3 replicates); genes with signal > 50 in at least half of the
    24 samples, log2."""
    x = pd.read_excel(os.path.join(DATA, 'kellermeier', 'Kellermeier2014_Supplemental_DS1.xlsx'), header=None)
    cols = [str(c) if pd.notna(c) else f'nan_{j}' for j, c in enumerate(x.iloc[2])]
    d = x.iloc[4:].copy(); d.columns = cols
    d = d[d['Media code'].notna()].copy()
    for c in d.columns[6:]:
        d[c] = pd.to_numeric(d[c], errors='coerce')
    media = d[['Media', 'nan_2', 'nan_3', 'nan_4']].astype(str).agg(''.join, axis=1)
    for f in 'NPKS':
        d[f] = media.str.contains(f).astype(int)
    d['L'] = (d['Day Length'].astype(str).str.startswith('Long')).astype(int)
    d.to_pickle(os.path.join(OUT, 'kellermeier_ds1.pkl'))
    print('Kellermeier DS1:', len(d), 'plants')

    x = pd.read_excel(os.path.join(DATA, 'kellermeier', 'Kellermeier2014_Supplemental_DS4.xlsx'),
                      sheet_name='RAW ABSOLUTE EXPRESSION', header=None)
    v = x.iloc[4:, 4:28].astype(float); v.index = x.iloc[4:, 2].astype(str).values
    keep = (v > 50).mean(axis=1) >= 0.5
    L = np.log2(v[keep].clip(lower=1)); L.columns = range(L.shape[1])
    D = pd.DataFrame({'N': (x.iloc[0, 4:28] == 'N').astype(int).values, 'P': (x.iloc[1, 4:28] == 'P').astype(int).values,
                      'K': (x.iloc[2, 4:28] == 'K').astype(int).values, 'rep': x.iloc[3, 4:28].values})
    L.to_pickle(os.path.join(OUT, 'kellermeier_ds4_log2.pkl')); D.to_pickle(os.path.join(OUT, 'kellermeier_ds4_design.pkl'))
    print('Kellermeier DS4:', L.shape[0], 'genes x', L.shape[1], 'samples')


# ----------------------------------------------------------------- Tekin et al. 2018 (E. coli, antibiotics)
CODE = re.compile(r'^(?:[A-Z]{3}\d){5}$')
NUM = re.compile(r'^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$')

def tekin():
    """The supplementary dataset (PDF) lists, for each five-antibiotic panel, a header '<code> <31 combination
    names>' followed by 3 replicate lines of 31 growth values (% of the no-drug control). The raw token stream
    is parsed; blocks broken by page layout are skipped."""
    raw = os.path.join(OUT, 'tekin_raw.txt')
    if not os.path.exists(raw):
        subprocess.run(['pdftotext', '-raw', os.path.join(DATA, 'tekin', 'Tekin2018_Dataset_MOESM2.pdf'), raw], check=True)
    toks = open(raw).read().split()
    blocks, i, n = [], 0, len(toks)
    while i < n:
        if not CODE.match(toks[i]):
            i += 1; continue
        code = toks[i]; names = toks[i + 1:i + 32]; drugs = re.findall(r'[A-Z]{3}', code)
        ok = len(names) == 31 and all(set(nm.split('+')) <= set(drugs) for nm in names)
        vals = toks[i + 32:i + 125]
        if ok and len(vals) == 93 and all(NUM.match(t) for t in vals):
            blocks.append((code, names, np.array(vals, float).reshape(3, 31))); i += 125
        else:
            i += 1
    pickle.dump(blocks, open(os.path.join(OUT, 'tekin_blocks.pkl'), 'wb'))
    print('Tekin:', len(blocks), 'panels')


def check(ref):
    for f in ['ristova_log2', 'ristova_design', 'kellermeier_ds1', 'kellermeier_ds4_log2', 'kellermeier_ds4_design']:
        a = pd.read_pickle(os.path.join(OUT, f + '.pkl')); b = pd.read_pickle(os.path.join(ref, f + '.pkl'))
        num = a.select_dtypes('number'); same = a.shape == b.shape and np.allclose(num.values, b[num.columns].values.astype(float), equal_nan=True)
        print(f, a.shape, 'identical to reference' if same else 'DIFFERS from reference')
    a = pickle.load(open(os.path.join(OUT, 'tekin_blocks.pkl'), 'rb')); b = pickle.load(open(os.path.join(ref, 'tekin_blocks.pkl'), 'rb'))
    same = len(a) == len(b) and all(x[0] == y[0] and np.allclose(x[2], y[2]) for x, y in zip(a, b))
    print('tekin_blocks', len(a), 'identical to reference' if same else 'DIFFERS from reference')


if __name__ == '__main__':
    ristova(); kellermeier(); tekin()
    if '--check' in sys.argv: check(sys.argv[sys.argv.index('--check') + 1])
