import numpy as np, pandas as pd, waves
waves.R = 1
GRID = [24, 48, 72, 96, 120, 144, 200, 264, 360, 480]
waves.CUTS = GRID
rng = np.random.default_rng(5)
A = np.load('fig/plan_K18_240.npy')                 # 240 first mixtures (same as before)
A = waves.augment(A, 480 - len(A), rng, cand=250)
np.save('fig/plan_K18_480.npy', A)
rows = []
for order, M, d, B in [(2, 5, 0.5, 0.1), (2, 5, 1.0, 0.1), (2, 5, 2.0, 0.1), (2, 5, 1.0, 0.45), (3, 3, 1.0, 0.1), (3, 3, 2.0, 0.1)]:
    r = waves.evaluate(A, order, M, d, B, reps=6 if order == 3 else 8); rows += r
    print(order, M, d, B, [(x['N'], round(x['recall'], 2), round(x['main'], 2)) for x in r], flush=True)
pd.DataFrame(rows).to_csv('fig/rep_sim_R1.csv', index=False)
