"""Figures and summary tables of the benchmark."""
import pickle, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
os.makedirs('fig', exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
COL = {'combinatorial': '#1F6F5C', 'random': '#E0A030', 'one_at_a_time': '#D1495B'}
LAB = {'combinatorial': 'Combinatoire (Shasha + foldover)', 'random': 'Aléatoire équilibré', 'one_at_a_time': 'Un à la fois (seuls, paires…)'}
DS = {'ristova': 'Ristova – transcriptome (15 348 gènes)', 'kellermeier': 'Kellermeier – traits racinaires (moitiés indépendantes)',
      'tekin': 'Tekin – croissance E. coli (2 000 blocs)'}

def load(name):
    p = f'res_{name}.pkl'
    return pd.DataFrame(pickle.load(open(p, 'rb'))) if os.path.exists(p) else None

def panel(ax, df, metric, title, ylabel, sim=True, ylim=(0, 1.02)):
    r = df[df.kind == 'real']
    for d, c in COL.items():
        g = r[r.design == d].groupby('N')[metric].mean()
        ax.plot(g.index, g.values, '-o', color=c, ms=3.5, lw=2, label=LAB[d])
        if sim and metric in df[df.kind == 'simulated'].columns:
            gs = df[(df.kind == 'simulated') & (df.design == d)].groupby('N')[metric].mean()
            ax.plot(gs.index, gs.values, '--', color=c, lw=1.2, alpha=.7)
    ax.set_title(title, fontsize=10, loc='left'); ax.set_xlabel('combinaisons mesurées (sur 32)'); ax.set_ylabel(ylabel)
    if ylim: ax.set_ylim(*ylim)
    ax.grid(alpha=.25); ax.set_xticks([6, 10, 16, 20, 24, 32])

if __name__ == '__main__':
    res = {k: load(k) for k in DS}
    res = {k: v for k, v in res.items() if v is not None}
    # Figure 1: main recall, interaction recall, false interactions, prediction error
    metrics = [('main_rec', 'Effets simples retrouvés', 'proportion', (0, 1.02)),
               ('int_rec', 'Interactions retrouvées', 'proportion', (0, 1.02)),
               ('false_int_eff', 'Réponses avec ≥ 1 fausse interaction', 'proportion', (0, 1.02)),
               ('pred_ratio', 'Erreur sur les combinaisons non mesurées', '× bruit d’une moyenne', None)]
    fig, axes = plt.subplots(len(res), 4, figsize=(15, 3.4 * len(res)), squeeze=False)
    for i, (k, df) in enumerate(res.items()):
        for j, (m, t, yl, lim) in enumerate(metrics):
            panel(axes[i, j], df, m, (DS[k] + '\n' if j == 0 else '\n') + t, yl, sim=(m != 'pred_ratio'), ylim=lim)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h + [plt.Line2D([], [], color='grey', ls='--')], l + ['simulation recalée sur les vrais effets'],
               loc='lower center', ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig('fig/benchmark_metrics.png', dpi=160, facecolor='white')
    # Figure 2: simulation vs reality (exact success)
    fig, ax = plt.subplots(figsize=(5.2, 5))
    mk = {'ristova': 'o', 'kellermeier': 's', 'tekin': '^'}
    for k, df in res.items():
        a = df[df.kind == 'real'].groupby(['design', 'N'])['success'].mean()
        b = df[df.kind == 'simulated'].groupby(['design', 'N'])['success'].mean()
        j = a.index.intersection(b.index)
        ax.scatter(b[j], a[j], marker=mk[k], s=28, color=[COL[d] for d, _ in j], label=DS[k].split(' –')[0], edgecolor='white', lw=.4)
    ax.plot([0, 1], [0, 1], color='grey', lw=1, ls='--'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel('réussite prédite par la simulation'); ax.set_ylabel('réussite observée sur données réelles')
    ax.set_title('Simulation recalée vs réalité\n(chaque point : un plan × un nombre de combinaisons)', fontsize=10, loc='left')
    ax.legend(frameon=False, fontsize=9); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig('fig/sim_vs_real.png', dpi=160, facecolor='white')
    # tables
    rows = []
    for k, df in res.items():
        r = df[df.kind == 'real']
        for d in COL:
            for N in (8, 10, 12, 16, 20, 24, 32):
                g = r[(r.design == d) & (r.N == N)]
                if len(g) == 0: continue
                rows.append(dict(dataset=k, design=d, N=N, success=g.success.mean(), agree_full=g.agree_full.mean(),
                                 main_rec=g.main_rec.mean(), int_rec=g.int_rec.mean(), false_int=g.false_int_eff.mean(),
                                 false_int_null=g.false_int_null.mean(), pred_ratio=g.pred_ratio.mean(),
                                 sim_success=df[(df.kind == 'simulated') & (df.design == d) & (df.N == N)].success.mean()))
    T = pd.DataFrame(rows); T.to_csv('fig/summary_table.csv', index=False)
    st = pd.DataFrame([dict(dataset=k, **df[df.kind == 'stop'][['stop_success', 'stop_N', 'stop_trap']].mean().to_dict(),
                            fixed32_success=df[(df.kind == 'real') & (df.N == 32)].success.mean()) for k, df in res.items()])
    st.to_csv('fig/stopping.csv', index=False)
    las = []
    for k, df in res.items():
        if 'lasso_int_exact' in df:
            g = df.dropna(subset=['lasso_int_exact']).groupby(['design', 'N'])[['fs_int_exact', 'lasso_int_exact']].mean()
            g['dataset'] = k; las.append(g.reset_index())
    if las: pd.concat(las).to_csv('fig/lasso.csv', index=False)
    pd.set_option('display.width', 220); print(T.round(2).to_string(index=False)); print(st.round(3)); print(pd.concat(las).round(2).to_string() if las else '')
