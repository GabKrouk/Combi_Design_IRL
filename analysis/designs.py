import itertools, numpy as np
from lib import all_cells
def covers(rows, t):
    for c in itertools.combinations(range(rows.shape[1]), t):
        if len({tuple(r) for r in rows[:, c]}) < 2 ** t: return False
    return True
def find_ca(K, t, N, seed=0, tries=200000):
    cells = all_cells(K); rng = np.random.default_rng(seed)
    for _ in range(tries):
        idx = rng.choice(len(cells), N, replace=False)
        if covers(cells[idx], t): return sorted(idx.tolist())
    return None
def ca_sizes(K):
    out = {}
    for t in (2, 3):
        for N in range(2 ** t, 2 ** K + 1):
            r = find_ca(K, t, N, tries=20000 if N < 2**K else 1)
            if r: out[t] = (N, r); break
    return out
def sequences(K, rng, ca3):
    cells = all_cells(K); ci = {tuple(c): i for i, c in enumerate(cells)}
    perm = rng.permutation(K); flip = rng.integers(0, 2, K)
    def tr(i):  # relabel factors + flip levels
        c = cells[i][perm] ^ flip; return ci[tuple(c)]
    ca = [tr(i) for i in ca3]
    fold = [ci[tuple(1 - cells[i])] for i in ca]
    rest = [i for i in rng.permutation(len(cells)) if i not in ca and i not in fold]
    combi = list(dict.fromkeys(ca + fold + rest))
    rand = list(rng.permutation(len(cells)))
    size = cells.sum(1); oat = []
    for s in range(K + 1): oat += list(rng.permutation(np.where(size == s)[0]))
    return dict(combinatorial=combi, random=rand, one_at_a_time=oat)
if __name__ == '__main__':
    print(ca_sizes(5))
