"""Figures for the manuscript (English, PNAS widths: 8.7 / 11.4 / 17.8 cm)."""
import os, pickle, numpy as np, pandas as pd, matplotlib
from math import comb
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
os.makedirs('paper', exist_ok=True)
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.3, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': .6, 'xtick.major.width': .6, 'ytick.major.width': .6,
                     'lines.linewidth': 1.4, 'savefig.dpi': 600, 'pdf.fonttype': 42})
CM = 1 / 2.54
TEAL, AMBER, RED, GREY, INK = '#1F6F5C', '#E0A030', '#D1495B', '#8C8C8C', '#222222'
DCOL = {'combinatorial': TEAL, 'random': AMBER, 'one_at_a_time': RED}
DLAB = {'combinatorial': 'Combinatorial design', 'random': 'Balanced random design', 'one_at_a_time': 'One-at-a-time ladder'}
DS = {'ristova': 'Arabidopsis transcriptome\n(5 signals, 15,348 genes)', 'kellermeier': 'Arabidopsis root traits\n(4 nutrients + light, 15 traits)',
      'tekin': 'E. coli growth\n(5 antibiotics, 2,000 panels)'}
DSHORT = {'Ristova': 'Plant transcriptome (5 signals)', 'Kellermeier_tx': 'Plant transcriptome (3 nutrients)',
          'Kellermeier_traits': 'Plant root traits', 'Tekin': 'Bacterial growth (antibiotics)'}
DSCOL = {'Ristova': TEAL, 'Kellermeier_tx': '#5FA8D3', 'Kellermeier_traits': AMBER, 'Tekin': RED}

def letter(ax, s, x=-0.18, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom', ha='left')

def save(fig, name):
    fig.savefig(f'paper/{name}.png', dpi=600, facecolor='white'); fig.savefig(f'paper/{name}.pdf', facecolor='white')
    plt.close(fig)

# ------------------------------------------------------------------ simulation summaries
def n_at(g, target):
    g = g.sort_values('N'); N = g.N.values.astype(float); r = np.maximum.accumulate(g.recall.values)
    if r[-1] < target: return np.nan
    i = int(np.argmax(r >= target))
    if i == 0: return N[0]
    f = (target - r[i - 1]) / max(r[i] - r[i - 1], 1e-9)
    return float(np.exp(np.log(N[i - 1]) + f * (np.log(N[i]) - np.log(N[i - 1]))))

def thresholds(df, keys):
    rows = []
    for k, g in df.groupby(keys):
        k = k if isinstance(k, tuple) else (k,)
        g95 = g.sort_values('N'); at = g95[g95.N >= (n_at(g, .95) or 0)].head(1)
        rows.append(dict(zip(keys, k), N80=n_at(g, .8), N95=n_at(g, .95),
                         fdr95=float(at.fdr.iloc[0]) if len(at) else np.nan,
                         main95=float(at.main_recall.iloc[0]) if len(at) else np.nan))
    return pd.DataFrame(rows)

def ex_pairs(K): return 1 + K + comb(K, 2)
def ex_trip(K): return 1 + K + comb(K, 2) + comb(K, 3)

P = pd.read_csv('fig/simK_pairs.csv'); TP = thresholds(P, ['K', 'M', 'd'])
T3 = pd.read_csv('fig/simK_trip.csv'); TT = thresholds(T3, ['K', 'M', 'd'])
BG = pd.read_csv('fig/simK_bg.csv') if os.path.exists('fig/simK_bg.csv') else None
TB = thresholds(BG, ['K', 'B', 'd']) if BG is not None and len(BG) else None
TP['ex'] = TP.K.map(ex_pairs); TP['gain'] = TP.ex / TP.N95
TT['ex'] = TT.K.map(ex_trip); TT['gain'] = TT.ex / TT.N95
TP.to_csv('paper/table_pairs_thresholds.csv', index=False); TT.to_csv('paper/table_triplets_thresholds.csv', index=False)
if TB is not None: TB.to_csv('paper/table_background_thresholds.csv', index=False)

# ================================================================== Figure 1: the idea and the cost
def fig1():
    fig = plt.figure(figsize=(17.8 * CM, 7.2 * CM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.55], wspace=.45, left=.05, right=.98, top=.86, bottom=.15)
    K = 10; rng = np.random.default_rng(3)
    # exhaustive pairs design
    ex = [np.zeros(K, int)] + [np.eye(K, dtype=int)[i] for i in range(K)]
    for i in range(K):
        for j in range(i + 1, K):
            v = np.zeros(K, int); v[[i, j]] = 1; ex.append(v)
    ex = np.array(ex)
    # combinatorial design: balanced, every pair in all four configurations
    best = None
    for _ in range(4000):
        A = np.zeros((16, K), int)
        for j in range(K): A[rng.permutation(16)[:8], j] = 1
        cnt = min(np.sum((A[:, i] == a) & (A[:, j] == b)) for i in range(K) for j in range(i + 1, K) for a in (0, 1) for b in (0, 1))
        if best is None or cnt > best[0]: best = (cnt, A)
    A = best[1]
    cmap = matplotlib.colors.ListedColormap(['#EDEDED', TEAL])
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1])
    for ax, M, ttl, n in ((axA, ex, 'Exhaustive pairwise screen', len(ex)), (axB, A, 'Combinatorial design', len(A))):
        ax.imshow(M, aspect='auto', cmap=cmap, interpolation='nearest', extent=(-.5, K - .5, n - .5, -.5))
        ax.set_ylim(len(ex) - .5, -.5)
        ax.set_xticks(range(K)); ax.set_xticklabels([str(i + 1) for i in range(K)], fontsize=5.5)
        ax.set_yticks([]); ax.set_xlabel('molecule'); ax.tick_params(length=0)
        for s in ax.spines.values(): s.set_visible(False)
        ax.set_title(f'{ttl}\n{n} conditions', fontsize=7)
        for x in np.arange(-.5, K, 1): ax.axvline(x, color='white', lw=.4)
        for y in np.arange(-.5, n, 1): ax.axhline(y, xmin=0, xmax=1, color='white', lw=.25)
    # highlight pair (2,6) in the combinatorial design
    i, j = 1, 5
    for r in range(len(A)):
        if A[r, i] and A[r, j]:
            for c in (i, j): axB.add_patch(Rectangle((c - .5, r - .5), 1, 1, fill=False, ec=AMBER, lw=1.6))
    for c in (i, j): axB.get_xticklabels()[c].set_color(AMBER); axB.get_xticklabels()[c].set_fontweight('bold')
    cnts = {(a, b): int(np.sum((A[:, i] == a) & (A[:, j] == b))) for a in (1, 0) for b in (1, 0)}
    axB.text(K / 2 - .5, len(A) + 2.2, 'molecules 2 & 6 together: %d×\n2 alone: %d×   6 alone: %d×   neither: %d×' %
             (cnts[(1, 1)], cnts[(1, 0)], cnts[(0, 1)], cnts[(0, 0)]), ha='center', va='top', fontsize=5.8, color=INK)
    axB.text(K / 2 - .5, len(A) + 9.5, 'every pair of molecules appears\nin all four presence/absence states', ha='center',
             va='top', fontsize=5.8, color=GREY, style='italic')
    letter(axA, 'A', -.08); letter(axB, 'B', -.08)
    # cost vs K
    ax = fig.add_subplot(gs[2]); Ks = np.arange(5, 101)
    ax.plot(Ks, [min(2.0 ** k, 1e31) for k in Ks], color=INK, lw=1.2, ls=':', label='Full factorial (all $2^K$ mixtures)')
    ax.plot(Ks, [ex_trip(k) for k in Ks], color=GREY, lw=1.4, ls='--', label='Exhaustive: all pairs and triplets')
    ax.plot(Ks, [ex_pairs(k) for k in Ks], color=GREY, lw=1.4, label='Exhaustive: all pairs')
    sub = TP[(TP.M == 5) & (TP.d == 1.0)].sort_values('K')
    ax.plot(sub.K, sub.N95, '-o', color=TEAL, ms=3.2, lw=1.8, label='Combinatorial design (this study)')
    sub3 = TT[(TT.M == 3) & (TT.d == 1.0)].sort_values('K')
    if sub3.N95.notna().any():
        ax.plot(sub3.K, sub3.N95, '-s', color=TEAL, ms=3, lw=1.2, mfc='white', label='   … when screening triplets too')
    for k in (18, 50, 100):
        r = sub[sub.K == k]
        if len(r) and np.isfinite(r.N95.iloc[0]):
            ax.annotate('%.0f× fewer\nmixtures' % (ex_pairs(k) / r.N95.iloc[0]), (k, r.N95.iloc[0]), xytext=(k - (4 if k == 100 else 0), r.N95.iloc[0] / 9),
                        ha='right' if k == 100 else 'center', fontsize=6, color=TEAL, fontweight='bold', arrowprops=dict(arrowstyle='-', color=TEAL, lw=.6))
    ax.set_yscale('log'); ax.set_ylim(8, 1e7); ax.set_xlim(4, 102)
    ax.set_xlabel('number of molecules tested (K)'); ax.set_ylabel('mixtures to measure')
    ax.set_title('Mixtures needed (5 interactions of size 1, 95% recall)', loc='left', fontsize=6.6)
    ax.legend(frameon=False, loc='upper left', handlelength=2.2); ax.grid(alpha=.2, which='major')
    letter(ax, 'C', -.16)
    save(fig, 'Fig1')

# ================================================================== Figure 2: scaling with K
def fig2():
    fig = plt.figure(figsize=(17.8 * CM, 11.8 * CM))
    gs = fig.add_gridspec(2, 3, wspace=.5, hspace=.62, left=.07, right=.99, top=.93, bottom=.1)
    axs = [fig.add_subplot(gs[0, i]) for i in range(3)] + [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1:])]
    ax = axs[0]; cols = {18: '#9BC9BC', 50: '#4E9C86', 100: '#123F35'}
    for K, c in cols.items():
        g = P[(P.K == K) & (P.M == 5) & (P.d == 1.0)].sort_values('N')
        ax.plot(g.N, g.recall, '-o', color=c, ms=2.6, label=f'K = {K}')
    ax.set_xscale('log'); ax.set_xlabel('mixtures measured (N)'); ax.set_ylabel('interactions recovered')
    ax.set_xticks([30, 100, 300]); ax.set_xticklabels(['30', '100', '300']); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_title('Recovering 5 pairwise interactions', loc='left'); ax.legend(frameon=False, loc='lower right'); ax.set_ylim(-.02, 1.03)
    ax.axhline(.95, color=GREY, lw=.5, ls=':'); letter(ax, 'A')
    ax = axs[1]; dc = {0.5: '#9BC9BC', 1.0: '#4E9C86', 2.0: '#123F35'}
    Ks = np.arange(10, 101); ax.plot(Ks, [ex_pairs(k) for k in Ks], color=GREY, label='exhaustive pairwise screen')
    for d, c in dc.items():
        s = TP[(TP.M == 5) & (TP.d == d)].sort_values('K'); ax.plot(s.K, s.N95, '-o', color=c, ms=2.6, label=f'interaction size {d:g}')
    ax.set_yscale('log'); ax.set_xlabel('molecules (K)'); ax.set_ylabel('mixtures for 95% recall')
    ax.set_title('Cost grows about linearly with K', loc='left'); ax.legend(frameon=False, loc='upper left', fontsize=5.8); letter(ax, 'B')
    ax = axs[2]; mc = {1: '#C7E3DA', 3: '#7FBCA9', 5: '#2F8A71', 10: '#123F35'}
    for M, c in mc.items():
        s = TP[(TP.M == M) & (TP.d == 1.0)].sort_values('K'); ax.plot(s.K, s.gain, '-o', color=c, ms=2.6, label=f'{M} interaction' + ('s' if M > 1 else ''))
    ax.set_yscale('log'); ax.set_xlabel('molecules (K)'); ax.set_ylabel('fold fewer mixtures\nthan exhaustive pairs')
    ax.set_title('Saving (interaction size 1)', loc='left'); ax.legend(frameon=False, loc='upper left', fontsize=5.8); letter(ax, 'C')
    ax.axhline(1, color=GREY, lw=.6)
    ax = axs[3]; mk = {1: 'o', 3: 's', 5: '^'}
    Ks = np.arange(10, 31); ax.plot(Ks, [ex_trip(k) for k in Ks], color=GREY, ls='--', label='exhaustive pairs + triplets')
    for M in (1, 3, 5):
        s = TT[(TT.M == M) & (TT.d == 1.0)].sort_values('K')
        ax.plot(s.K, s.N95, '-' + mk[M], color='#2F8A71', ms=2.8, lw=1, label=f'{M} three-way interaction' + ('s' if M > 1 else ''))
    ax.set_yscale('log'); ax.set_xlabel('molecules (K)'); ax.set_ylabel('mixtures for 95% recall')
    ax.set_title('Three-way interactions (size 1)', loc='left'); ax.legend(frameon=False, loc='upper left', fontsize=5.6); letter(ax, 'D')
    ax = axs[4]
    EX = pd.read_csv('fig/simK_exhaustive.csv'); x = 0; ticks = []; tl = []
    for K in (18, 30, 50):
        for d in (0.5, 1.0, 2.0):
            e = EX[(EX.K == K) & (EX.d == d)]
            c = TP[(TP.K == K) & (TP.M == 5) & (TP.d == d)]
            if not len(e) or not len(c): continue
            ax.bar(x - .2, e.recall.iloc[0], .38, color=GREY)
            ax.bar(x + .2, .95, .38, color=TEAL)
            ax.text(x - .2, e.recall.iloc[0] + .03, '%d' % e.N.iloc[0], ha='center', fontsize=5.2, color=INK, rotation=90, va='bottom')
            ax.text(x + .2, .98, '%d' % round(c.N95.iloc[0]), ha='center', fontsize=5.2, color=TEAL, rotation=90, va='bottom')
            ticks.append(x); tl.append(f'size {d:g}'); x += 1
        ax.text(x - 2, -.36, f'K = {K} molecules', ha='center', fontsize=6.3, transform=ax.transData)
        x += .6
    ax.set_xticks(ticks); ax.set_xticklabels(tl, fontsize=5.8); ax.set_ylim(0, 1.32); ax.set_yticks([0, .25, .5, .75, 1])
    ax.set_ylabel('interactions recovered')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GREY, label='exhaustive pairwise screen'), Patch(color=TEAL, label='combinatorial design (at its 95% point)')],
              frameon=False, loc='upper left', ncol=2, fontsize=6)
    ax.set_title('Same replication (3 per mixture): numbers above bars = mixtures measured', loc='left', fontsize=6.8); letter(ax, 'E', -.08)
    save(fig, 'Fig2')

# ================================================================== Figure 3: real full factorials
def load_res(n):
    return pd.DataFrame(pickle.load(open(f'res_{n}.pkl', 'rb')))

def fig3():
    res = {k: load_res(k) for k in DS}
    fig, axs = plt.subplots(3, 3, figsize=(17.8 * CM, 13.5 * CM), gridspec_kw=dict(hspace=.55, wspace=.35, left=.08, right=.99, top=.9, bottom=.1))
    rows = [('main_rec', 'single-molecule effects\nrecovered', (0, 1.03)), ('pred_ratio', 'error on untested mixtures\n(× replicate noise)', None),
            ('int_rec', 'interactions recovered', (0, 1.03))]
    for j, (k, df) in enumerate(res.items()):
        r = df[df.kind == 'real']
        for i, (m, yl, lim) in enumerate(rows):
            ax = axs[i, j]
            for d, c in DCOL.items():
                g = r[r.design == d].groupby('N')[m].mean()
                ax.plot(g.index, g.values, '-o', color=c, ms=2.4, label=DLAB[d])
            if lim: ax.set_ylim(*lim)
            if m == 'pred_ratio': ax.axhline(1, color=GREY, lw=.6, ls=':'); ax.set_ylim(0.8, 4)
            ax.set_xticks([6, 10, 16, 24, 32]); ax.grid(alpha=.2)
            if j == 0: ax.set_ylabel(yl)
            if i == 2: ax.set_xlabel('mixtures measured (of 32)')
            if i == 0:
                ax.set_title(DS[k], fontsize=7)
                sec = ax.secondary_xaxis('top', functions=(lambda x: x / 32 * 100, lambda x: x * 32 / 100))
                sec.set_xticks([25, 50, 100]); sec.set_xticklabels(['25%', '50%', '100%'], fontsize=5.6); sec.tick_params(length=2)
    for i, L in enumerate('ABC'): letter(axs[i, 0], L, -.32)
    h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc='lower center', ncol=3, frameon=False, fontsize=6.5)
    save(fig, 'Fig3')

# ================================================================== Figure 4: rarity of high-order interactions
def fig4():
    R = pd.read_csv('fig/rarity.csv'); RB = pd.read_csv('fig/recall_by_order.csv')
    order_dsets = ['Ristova', 'Kellermeier_tx', 'Kellermeier_traits', 'Tekin']
    fig, axs = plt.subplots(2, 2, figsize=(17.8 * CM, 11.5 * CM), gridspec_kw=dict(wspace=.32, hspace=.62, left=.08, right=.98, top=.93, bottom=.13))
    ax = axs[0, 0]; oc = {1: TEAL, 2: AMBER, 3: RED, 4: '#7B3F61', 5: '#3B1F2B'}
    on = {1: 'single molecules', 2: 'pairs', 3: 'triplets', 4: 'four-way', 5: 'five-way'}
    for x, dname in enumerate(order_dsets):
        bot = 0
        for o in range(1, 6):
            v = R[(R.dataset == dname) & (R.order == o)].var_share_pooled
            if len(v):
                ax.bar(x, v.iloc[0] * 100, bottom=bot, color=oc[o], width=.65, label=on[o] if x == 0 else None, edgecolor='white', lw=.3)
                if v.iloc[0] > .04: ax.text(x, bot + v.iloc[0] * 50, '%.0f%%' % (v.iloc[0] * 100), ha='center', va='center', fontsize=5.6, color='white')
                bot += v.iloc[0] * 100
    ax.set_xticks(range(4)); ax.set_xticklabels(['Transcriptome\n(5 signals)', 'Transcriptome\n(3 nutrients)', 'Root traits', 'Bacterial\ngrowth'], fontsize=5.8); ax.set_xlim(-.5, 5.3)
    ax.set_ylabel('% of variation among mixtures'); ax.set_title('Share of the response explained by each order', loc='left')
    ax.legend(frameon=False, fontsize=6, loc='center right', handlelength=1); letter(ax, 'A', -.13)
    ax = axs[0, 1]
    for dname in order_dsets:
        g = R[R.dataset == dname].sort_values('order'); per = g.var_share_pooled / g.n_terms
        ax.plot(g.order, per * 100, '-o', color=DSCOL[dname], ms=3, label=DSHORT[dname])
    ax.set_yscale('log'); ax.set_xticks([1, 2, 3, 4, 5]); ax.set_xlabel('order of the effect'); ax.set_ylabel('% of variation per term')
    ax.set_title('Weight of an individual term', loc='left'); letter(ax, 'B', -.13)
    h, l = ax.get_legend_handles_labels(); fig.legend(h, l, loc='lower center', ncol=4, frameon=False, fontsize=6.3)
    ax = axs[1, 0]
    for dname in order_dsets:
        g = R[R.dataset == dname].sort_values('order'); ax.plot(g.order, g.frac_resp_with_sig * 100, '-o', color=DSCOL[dname], ms=3)
    ax.set_xticks([1, 2, 3, 4, 5]); ax.set_xlabel('order of the effect'); ax.set_ylabel('% of responsive outputs\nwith ≥ 1 detectable term')
    ax.set_title('How often each order is detectable (full design)', loc='left'); letter(ax, 'C', -.13)
    ax = axs[1, 1]; rm = {'Ristova': 'Ristova', 'Kellermeier': 'Kellermeier_traits', 'Tekin': 'Tekin'}
    r = RB[(RB.plan == 'combinatorial') & (RB.type == 'réel')].groupby(['jeu', 'ordre', 'N']).rappel.mean().reset_index()
    for jeu, dname in rm.items():
        for o, ls, m in ((2, '-', 'o'), (3, '--', 's')):
            g = r[(r.jeu == jeu) & (r.ordre == o)]; ax.plot(g.N, g.rappel, ls + m, color=DSCOL[dname], ms=3, lw=1.2)
    ax.plot([], [], '-o', color=GREY, ms=3, label='pairs'); ax.plot([], [], '--s', color=GREY, ms=3, label='triplets')
    ax.set_xticks([10, 16, 24, 32]); ax.set_xlabel('mixtures measured (of 32)'); ax.set_ylabel('interactions recovered')
    ax.set_title('Triplets are harder to recover than pairs', loc='left'); ax.set_ylim(-.02, 1.03); letter(ax, 'D', -.13)
    ax.legend(frameon=False, fontsize=6, loc='upper left')
    save(fig, 'Fig4')

# ================================================================== Figure 5: calibration
def fig5():
    C = pd.read_csv('fig/recalibration_check.csv')
    C = C[~((C.jeu == 'Ristova') & (C.ordre == 3))]
    fig, axs = plt.subplots(1, 3, figsize=(17.8 * CM, 6.2 * CM), gridspec_kw=dict(wspace=.42, left=.07, right=.99, top=.84, bottom=.2))
    mk = {'Ristova': 'o', 'Kellermeier traits': 's', 'Tekin': '^'}; lab = {'Ristova': 'transcriptome', 'Kellermeier traits': 'root traits', 'Tekin': 'bacterial growth'}
    for ax, col, ttl, L in ((axs[0], 'ancien_sim', 'Naive simulation\n(true effects + Gaussian noise)', 'A'),
                            (axs[1], 'recalé', 'Calibrated simulation\n(empirical effect-size model)', 'B')):
        for j, g in C.groupby('jeu'):
            ax.scatter(g[col], g['réel'], marker=mk[j], s=12, c=[TEAL if o == 2 else RED for o in g.ordre], edgecolor='white', lw=.3)
        ax.plot([0, 1], [0, 1], color=GREY, lw=.7, ls='--'); ax.set_xlim(-.02, 1.02); ax.set_ylim(-.02, 1.02)
        b = (C[col] - C['réel']).mean(); a = (C[col] - C['réel']).abs().mean()
        ax.text(.04, .96, f'mean bias {b:+.2f}\nmean |error| {a:.2f}', transform=ax.transAxes, va='top', fontsize=6)
        ax.set_xlabel('recall predicted by simulation'); ax.set_ylabel('recall observed in real data'); ax.set_title(ttl, loc='left', fontsize=6.8)
        letter(ax, L, -.25)
    for j, m in mk.items(): axs[0].scatter([], [], marker=m, c='white', edgecolor=INK, lw=.6, s=14, label=lab[j])
    axs[0].scatter([], [], c=TEAL, s=12, label='pairs'); axs[0].scatter([], [], c=RED, s=12, label='triplets')
    axs[0].legend(frameon=False, fontsize=5.6, loc='lower right')
    ax = axs[2]
    if TB is not None and len(TB):
        for K, c in ((18, '#9BC9BC'), (50, '#4E9C86'), (100, '#123F35')):
            for B, ls, mkk in ((0.0, '-', 'o'), (0.45, '--', 's')):
                s = TB[(TB.K == K) & (TB.B == B)].sort_values('d')
                ax.plot(s.d, s.N95, ls + mkk, color=c, ms=2.6, lw=1.1, label=f'K={K}' if B == 0 else None)
        ax.plot([], [], '-', color=GREY, label='sparse (no background)'); ax.plot([], [], '--', color=GREY, label='dense small interactions')
        ax.set_yscale('log'); ax.set_xscale('log'); ax.set_xticks([.5, 1, 2]); ax.set_xticklabels(['0.5', '1', '2']); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_yticks([60, 100, 200, 400]); ax.set_yticklabels(['60', '100', '200', '400'])
        ax.set_xlabel('size of the interactions sought\n(× replicate SD)'); ax.set_ylabel('N for 95% recall (5 interactions)')
        ax.set_title('What drives the cost', loc='left'); ax.legend(frameon=False, fontsize=5.4, loc='upper right')
    letter(ax, 'C', -.25)
    save(fig, 'Fig5')

# ================================================================== SI figures
def figS1():
    rows = []
    for k in DS:
        df = load_res(k); st = df[df.kind == 'stop']
        rows.append(dict(ds=k, succ=st.stop_success.mean(), N=st.stop_N.mean(), trap=st.stop_trap.mean(),
                         full=df[(df.kind == 'real') & (df.N == 32) & (df.design == 'combinatorial')].success.mean()))
    R = pd.DataFrame(rows)
    fig, axs = plt.subplots(1, 3, figsize=(17.8 * CM, 5.5 * CM), gridspec_kw=dict(wspace=.45, left=.07, right=.99, top=.85, bottom=.2))
    x = np.arange(3); lab = ['Plant\ntranscriptome', 'Plant\nroot traits', 'Bacterial\ngrowth']
    ax = axs[0]; ax.bar(x - .2, R.full, .38, color=GREY, label='all 32 mixtures'); ax.bar(x + .2, R.succ, .38, color=TEAL, label='wave rule')
    ax.set_xticks(x); ax.set_xticklabels(lab, rotation=0, ha='center'); ax.set_ylabel('outputs fully recovered'); ax.legend(frameon=False)
    ax.set_title('Conclusions reached', loc='left'); letter(ax, 'A')
    ax = axs[1]; ax.bar(x, R.N, .5, color=TEAL); ax.axhline(32, color=GREY, ls=':', lw=.8)
    ax.set_xticks(x); ax.set_xticklabels(lab, rotation=0, ha='center'); ax.set_ylabel('mixtures used (mean)'); ax.set_title('Cost', loc='left'); letter(ax, 'B')
    ax = axs[2]; ax.bar(x, R.trap * 100, .5, color=RED); ax.set_xticks(x); ax.set_xticklabels(lab, rotation=0, ha='center')
    ax.set_ylabel('% stopped early and wrong'); ax.set_title('Premature stops', loc='left'); letter(ax, 'C')
    save(fig, 'FigS1'); R.to_csv('paper/table_stopping.csv', index=False)

def figS2():
    L = pd.read_csv('fig/lasso.csv')
    fig, axs = plt.subplots(1, 2, figsize=(11.4 * CM, 5.2 * CM), gridspec_kw=dict(wspace=.35, left=.1, right=.99, top=.85, bottom=.2))
    for ax, k, L_ in zip(axs, ('ristova', 'tekin'), 'AB'):
        g = L[(L.dataset == k) & (L.design == 'combinatorial')]
        ax.plot(g.N, g.fs_int_exact, '-o', color=TEAL, ms=2.6, label='forward selection (this study)')
        ax.plot(g.N, g.lasso_int_exact, '-s', color=AMBER, ms=2.6, label='lasso (BIC)')
        ax.set_ylim(0, 1.03); ax.set_xlabel('mixtures measured (of 32)'); ax.set_ylabel('exact set of interactions')
        ax.set_title(DS[k].split('\n')[0], loc='left'); letter(ax, L_, -.22)
    axs[0].legend(frameon=False, loc='lower left')
    save(fig, 'FigS2')

def figS3():
    fig, axs = plt.subplots(2, 3, figsize=(17.8 * CM, 10 * CM), sharey=True, gridspec_kw=dict(wspace=.12, hspace=.45, left=.07, right=.99, top=.9, bottom=.1))
    dc = {0.5: '#9BC9BC', 1.0: '#4E9C86', 2.0: '#123F35'}; ls = {1: ':', 3: '-.', 5: '-', 10: '--'}
    for ax, K in zip(axs.ravel(), (10, 18, 30, 50, 75, 100)):
        for (M, d), g in P[P.K == K].groupby(['M', 'd']):
            g = g.sort_values('N'); ax.plot(g.N, g.recall, ls[M], color=dc[d], lw=1)
        ax.axvline(ex_pairs(K), color=GREY, lw=.7); ax.set_xscale('log'); ax.set_title(f'K = {K} molecules', loc='left')
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xticks([20, 50, 100, 200, 500, 1000, 5000]); ax.set_xticklabels(['20', '50', '100', '200', '500', '1000', '5000'], fontsize=5.5)
        ax.set_xlim(min(P[P.K == K].N) * .9, ex_pairs(K) * 1.3)
    for a in axs[:, 0]: a.set_ylabel('interactions recovered')
    for a in axs[1]: a.set_xlabel('mixtures measured (N)')
    h = [plt.Line2D([], [], color=c, lw=1.4, label=f'size {d:g}') for d, c in dc.items()] + \
        [plt.Line2D([], [], color=GREY, ls=l, lw=1, label=f'{M} interaction' + ('s' if M > 1 else '')) for M, l in ls.items()] + \
        [plt.Line2D([], [], color=GREY, lw=.7, label='exhaustive pairs')]
    fig.legend(handles=h, loc='upper center', ncol=8, frameon=False, fontsize=5.8)
    save(fig, 'FigS3')


def figS4():
    S3 = pd.read_csv('fig/rep_sim_R3.csv'); S1 = pd.read_csv('fig/rep_sim_R1.csv')
    fig, axs = plt.subplots(1, 2, figsize=(11.4 * CM, 5.2 * CM), gridspec_kw=dict(wspace=.35, left=.11, right=.99, top=.86, bottom=.2))
    for ax, (order, tit, L) in zip(axs, ((2, 'Five pairwise interactions', 'A'), (3, 'Three three-way interactions', 'B'))):
        a = S3[(S3.order == order) & (S3.d == 1) & (S3.B == .1)]; b = S1[(S1.order == order) & (S1.d == 1) & (S1.B == .1)]
        ax.plot(a.N * 3, a.recall, '-s', color=GREY, ms=2.6, label='3 replicates per mixture')
        ax.plot(b.N, b.recall, '-o', color=TEAL, ms=2.6, label='1 measurement per mixture')
        ax.set_xlim(0, 500); ax.set_ylim(-.02, 1.05); ax.set_title(tit + ' (size 1)', loc='left'); ax.grid(alpha=.2)
        ax.set_xlabel('total measurements'); ax.set_ylabel('interactions recovered'); letter(ax, L, -.25)
    axs[0].legend(frameon=False, loc='lower right', fontsize=5.8)
    save(fig, 'FigS4')

def figS5():
    D = pd.read_csv('fig/simK_noise.csv'); C = pd.read_csv('fig/slope_by_noise.csv'); R = pd.read_csv('fig/rule_ci.csv')
    fig, axs = plt.subplots(1, 2, figsize=(11.4 * CM, 5.4 * CM), gridspec_kw=dict(wspace=.45, left=.11, right=.99, top=.86, bottom=.2))
    col = {'gaussian': GREY, 'ristova': TEAL, 'kellermeier': AMBER, 'tekin': RED}
    lab = {'gaussian': 'Gaussian noise', 'ristova': 'noise of the transcriptome', 'kellermeier': 'noise of the root traits', 'tekin': 'noise of bacterial growth'}
    ax = axs[0]
    for nz, g in D.groupby('noise'):
        pts = [(K, n_at(h.rename(columns={}), .95)) for K, h in g.assign(fdr=0, main_recall=1).groupby('K')]
        pts = np.array(pts); ax.plot(pts[:, 0], pts[:, 1], '-o', color=col[nz], ms=2.6, label=lab[nz])
    ax.set_xlabel('molecules (K)'); ax.set_ylabel('mixtures for 95% recall'); ax.set_title('Empirical noise', loc='left')
    ax.legend(frameon=False, fontsize=5.4, loc='upper left'); letter(ax, 'A', -.25)
    ax = axs[1]; rows = []
    k = R[(R.term == 'K') & (R.method == 'scenario + Monte Carlo')].iloc[0]; rows.append(('all scenarios', k.estimate, k.lo, k.hi, INK))
    for nz in ('gaussian', 'ristova', 'kellermeier', 'tekin'):
        c = C[C.noise == nz].iloc[0]; rows.append((lab[nz].replace('noise of the ', '').replace('noise of ', ''), c.slope, c.lo, c.hi, col[nz]))
    for i, (n, e, lo, hi, c) in enumerate(rows):
        ax.errorbar(e, i, xerr=[[e - lo], [hi - e]], fmt='o', color=c, ms=3.5, capsize=2, lw=1)
    ax.axvline(1.618, color=GREY, ls=':', lw=.7); ax.text(1.625, len(rows) - .6, 'φ', fontsize=7, color=GREY)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=6); ax.invert_yaxis()
    ax.set_xlabel('additional mixtures per molecule'); ax.set_title('Coefficient of K (95% CI)', loc='left'); letter(ax, 'B', -.55)
    save(fig, 'FigS5')

if __name__ == '__main__':
    import sys
    for f in (sys.argv[1:] or ['fig1', 'fig2', 'fig3', 'fig4', 'fig5']): globals()[f]()
    pd.set_option('display.width', 200)
    print(TP.round(2).to_string(index=False)); print(TT.round(2).to_string(index=False))
    if TB is not None: print(TB.round(2).to_string(index=False))

