"""Share of the variance of the true combination means carried by each interaction order (noise-corrected),
and fraction of responses with >=1 significant term of each order (q<0.01)."""
import pickle, numpy as np, pandas as pd
out, dists = pickle.load(open('prep/truth.pkl', 'rb'))
rows = []
for name, tr in out.items():
    order = np.array([len(t) for t in tr['terms']]); K = order.max()
    sigq = tr['q'] < 0.01                                   # G x T
    eff = sigq.any(1)
    c2 = np.clip(tr['coef'] ** 2 - tr['se'] ** 2, 0, None).T  # G x T, noise-corrected
    tot = c2.sum(1)
    for o in range(1, K + 1):
        m = order == o
        share = (c2[:, m].sum(1) / np.where(tot > 0, tot, 1))[eff]
        rows.append(dict(dataset=name, order=o, n_terms=int(m.sum()),
                         var_share=float(np.mean(share)), var_share_pooled=float(c2[eff][:, m].sum() / tot[eff].sum()),
                         frac_resp_with_sig=float(sigq[eff][:, m].any(1).mean()),
                         frac_terms_sig=float(sigq[eff][:, m].mean()),
                         median_abs_sig=float(np.median(np.abs(tr['coef'].T[sigq & m[None, :]] / tr['sigma'][:, None].repeat(m.size, 1)[sigq & m[None, :]]) * 2)) if (sigq & m[None, :]).any() else np.nan,
                         n_resp=int(eff.sum())))
df = pd.DataFrame(rows); df.to_csv('fig/rarity.csv', index=False)
pd.set_option('display.width', 200); print(df.round(4).to_string(index=False))
