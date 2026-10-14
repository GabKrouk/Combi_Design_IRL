"""Fill numbers, resolve citations (PNAS numeric style), build a JSON block list for the docx builder."""
import re, json, numpy as np, pandas as pd
from math import comb

REFS = {
 'ristova2016': 'D. Ristova *et al.*, Combinatorial interaction network of transcriptomic and phenotypic responses to nitrogen and hormones in the *Arabidopsis thaliana* root. *Sci. Signal.* **9**, rs13 (2016).',
 'kellermeier2014': 'F. Kellermeier *et al.*, Analysis of the root system architecture of *Arabidopsis* provides a quantitative readout of crosstalk between nutritional signals. *Plant Cell* **26**, 1480–1496 (2014).',
 'yeh2006': 'P. Yeh, A. I. Tschumi, R. Kishony, Functional classification of drugs by properties of their pairwise interactions. *Nat. Genet.* **38**, 489–494 (2006).',
 'beppler2016': 'C. Beppler *et al.*, Uncovering emergent interactions in three-way combinations of stressors. *J. R. Soc. Interface* **13**, 20160800 (2016).',
 'tekin2018': 'E. Tekin *et al.*, Prevalence and patterns of higher-order drug interactions in *Escherichia coli*. *NPJ Syst. Biol. Appl.* **4**, 31 (2018).',
 'lehar2009': 'J. Lehár *et al.*, Synergistic drug combinations tend to improve therapeutically relevant selectivity. *Nat. Biotechnol.* **27**, 659–666 (2009).',
 'holbeck2017': 'S. L. Holbeck *et al.*, The National Cancer Institute ALMANAC: A comprehensive screening resource for the detection of anticancer drug pairs with enhanced therapeutic activity. *Cancer Res.* **77**, 3564–3576 (2017).',
 'friedman2017': 'J. Friedman, L. M. Higgins, J. Gore, Community structure follows simple assembly rules in microbial microcosms. *Nat. Ecol. Evol.* **1**, 0109 (2017).',
 'skwara2023': 'A. Skwara *et al.*, Statistically learning the functional landscape of microbial communities. *Nat. Ecol. Evol.* **7**, 1823–1833 (2023).',
 'bliss1939': 'C. I. Bliss, The toxicity of poisons applied jointly. *Ann. Appl. Biol.* **26**, 585–615 (1939).',
 'fisher1935': 'R. A. Fisher, *The Design of Experiments* (Oliver and Boyd, Edinburgh, 1935).',
 'plackett1946': 'R. L. Plackett, J. P. Burman, The design of optimum multifactorial experiments. *Biometrika* **33**, 305–325 (1946).',
 'box2005': 'G. E. P. Box, J. S. Hunter, W. G. Hunter, *Statistics for Experimenters: Design, Innovation, and Discovery* (Wiley, Hoboken, NJ, ed. 2, 2005).',
 'wu2009': 'C. F. J. Wu, M. S. Hamada, *Experiments: Planning, Analysis, and Optimization* (Wiley, Hoboken, NJ, ed. 2, 2009).',
 'hamada1992': 'M. Hamada, C. F. J. Wu, Analysis of designed experiments with complex aliasing. *J. Qual. Technol.* **24**, 130–137 (1992).',
 'chipman1996': 'H. Chipman, Bayesian variable selection with related predictors. *Can. J. Stat.* **24**, 17–36 (1996).',
 'li2006': 'X. Li, N. Sudarsanam, D. D. Frey, Regularities in data from factorial experiments. *Complexity* **11**, 32–45 (2006).',
 'sloane1993': 'N. J. A. Sloane, Covering arrays and intersecting codes. *J. Comb. Des.* **1**, 51–63 (1993).',
 'cohen1997': 'D. M. Cohen, S. R. Dalal, M. L. Fredman, G. C. Patton, The AETG system: An approach to testing based on combinatorial design. *IEEE Trans. Softw. Eng.* **23**, 437–444 (1997).',
 'colbourn2004': 'C. J. Colbourn, Combinatorial aspects of covering arrays. *Le Matematiche* **59**, 125–172 (2004).',
 'shasha2001': 'D. E. Shasha, A. Y. Kouranov, L. V. Lejay, M. F. Chou, G. M. Coruzzi, Using combinatorial design to study regulation by multiple input signals. A tool for parsimony in the post-genomics era. *Plant Physiol.* **127**, 1590–1594 (2001).',
 'dorfman1943': 'R. Dorfman, The detection of defective members of large populations. *Ann. Math. Stat.* **14**, 436–440 (1943).',
 'du2000': 'D.-Z. Du, F. K. Hwang, *Combinatorial Group Testing and Its Applications* (World Scientific, Singapore, ed. 2, 2000).',
 'angluin2008': 'D. Angluin, J. Chen, Learning a hidden graph using O(log n) queries per edge. *J. Comput. Syst. Sci.* **74**, 546–556 (2008).',
 'colbourn2008': 'C. J. Colbourn, D. W. McClary, Locating and detecting arrays for interaction faults. *J. Comb. Optim.* **15**, 17–48 (2008).',
 'donoho2006': 'D. L. Donoho, Compressed sensing. *IEEE Trans. Inf. Theory* **52**, 1289–1306 (2006).',
 'poelwijk2019': 'F. J. Poelwijk, M. Socolich, R. Ranganathan, Learning the pattern of epistasis linking genotype and phenotype in a protein. *Nat. Commun.* **10**, 4213 (2019).',
 'brookes2022': 'D. H. Brookes, A. Aghazadeh, J. Listgarten, On the sparsity of fitness functions and implications for learning. *Proc. Natl. Acad. Sci. U.S.A.* **119**, e2109649118 (2022).',
 'wood2012': 'K. Wood, S. Nishida, E. D. Sontag, P. Cluzel, Mechanism-independent method for predicting response to multidrug combinations in bacteria. *Proc. Natl. Acad. Sci. U.S.A.* **109**, 12254–12259 (2012).',
 'zimmer2016': 'A. Zimmer, I. Katzir, E. Dekel, A. E. Mayo, U. Alon, Prediction of multidimensional drug dose responses based on measurements of drug pairs. *Proc. Natl. Acad. Sci. U.S.A.* **113**, 10442–10447 (2016).',
 'zimmer2017': 'A. Zimmer, A. Tendler, I. Katzir, A. Mayo, U. Alon, Prediction of drug cocktail effects when the number of measurements is limited. *PLoS Biol.* **15**, e2002518 (2017).',
 'katzir2019': 'I. Katzir, M. Cokol, B. B. Aldridge, U. Alon, Prediction of ultra-high-order antibiotic combinations based on pairwise interactions. *PLoS Comput. Biol.* **15**, e1006774 (2019).',
 'russ2018': 'D. Russ, R. Kishony, Additivity of inhibitory effects in multidrug combinations. *Nat. Microbiol.* **3**, 1339–1345 (2018).',
 'kleitman1973': 'D. J. Kleitman, J. Spencer, Families of k-independent sets. *Discrete Math.* **6**, 255–262 (1973).',
 'tibshirani1996': 'R. Tibshirani, Regression shrinkage and selection via the lasso. *J. R. Stat. Soc. B* **58**, 267–288 (1996).',
 'efron2004': 'B. Efron, T. Hastie, I. Johnstone, R. Tibshirani, Least angle regression. *Ann. Stat.* **32**, 407–499 (2004).',
 'weinreich2013': 'D. M. Weinreich, Y. Lan, C. S. Wylie, R. B. Heckendorn, Should evolutionary geneticists worry about higher-order epistasis? *Curr. Opin. Genet. Dev.* **23**, 700–707 (2013).',
 'sailer2017': 'Z. R. Sailer, M. J. Harms, Detecting high-order epistasis in nonlinear genotype-phenotype maps. *Genetics* **205**, 1079–1088 (2017).',
 'otwinowski2018': 'J. Otwinowski, D. M. McCandlish, J. B. Plotkin, Inferring the shape of global epistasis. *Proc. Natl. Acad. Sci. U.S.A.* **115**, E7550–E7558 (2018).',
 'bien2013': 'J. Bien, J. Taylor, R. Tibshirani, A lasso for hierarchical interactions. *Ann. Stat.* **41**, 1111–1141 (2013).',
 'lim2015': 'M. Lim, T. Hastie, Learning interactions via hierarchical group-lasso regularization. *J. Comput. Graph. Stat.* **24**, 627–654 (2015).',
 'gelman2007': 'A. Gelman, J. Hill, *Data Analysis Using Regression and Multilevel/Hierarchical Models* (Cambridge University Press, Cambridge, 2007).',
 'stephens2017': 'M. Stephens, False discovery rates: A new deal. *Biostatistics* **18**, 275–294 (2017).',
 'edgar2002': 'R. Edgar, M. Domrachev, A. E. Lash, Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Res.* **30**, 207–210 (2002).',
 'benjamini1995': 'Y. Benjamini, Y. Hochberg, Controlling the false discovery rate: A practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995).',
 'efron2010': 'B. Efron, *Large-Scale Inference: Empirical Bayes Methods for Estimation, Testing, and Prediction* (Cambridge University Press, Cambridge, 2010).',
 'pedregosa2011': 'F. Pedregosa *et al.*, Scikit-learn: Machine learning in Python. *J. Mach. Learn. Res.* **12**, 2825–2830 (2011).',
}

B = '../analysis/'
TP = pd.read_csv(B + 'paper/table_pairs_thresholds.csv'); TB = pd.read_csv(B + 'paper/table_background_thresholds.csv')
EX = pd.read_csv(B + 'fig/simK_exhaustive.csv')
def tp(K, M, d, c='N95'): return float(TP[(TP.K == K) & (TP.M == M) & (TP.d == d)][c].iloc[0])
V = dict(
    main95_min='%.0f' % (100 * TP.main95.min()), fdr95_mean='%.0f' % (100 * TP.fdr95.mean()), fdr95_max='%.0f' % (100 * TP.fdr95.max()),
    N_K18='%.0f' % tp(18, 5, 1.0), N_K50='%.0f' % tp(50, 5, 1.0), N_K100='%.0f' % tp(100, 5, 1.0),
    gain18='%.0f' % tp(18, 5, 1.0, 'gain'), gain50='%.0f' % tp(50, 5, 1.0, 'gain'), gain100='%.0f' % tp(100, 5, 1.0, 'gain'),
    N_K100_d05='%.0f' % tp(100, 5, 0.5), N_K100_M1='%.0f' % tp(100, 1, 1.0), N_K100_M10='%.0f' % tp(100, 10, 1.0),
    ex_d1_min='%.0f' % (100 * EX[(EX.d == 1) & (EX.K >= 18)].recall.min()), ex_d1_max='%.0f' % (100 * EX[(EX.d == 1) & (EX.K >= 18)].recall.max()),
    bg_N50_B0='%.0f' % float(TB[(TB.K == 50) & (TB.B == 0) & (TB.d == 1)].N95.iloc[0]),
    bg_N50_B45='%.0f' % float(TB[(TB.K == 50) & (TB.B == .45) & (TB.d == 1)].N95.iloc[0]),
)

def parse(path):
    """Very small markdown: % meta, # h1, ## h2, paragraphs."""
    meta, blocks = {}, []
    for chunk in open(path).read().split('\n\n'):
        for line in [l for l in chunk.split('\n') if l.startswith('%')]:
            k, v = line[1:].split(' ', 1); meta[k] = v
        body = '\n'.join(l for l in chunk.split('\n') if not l.startswith('%')).strip()
        if not body: continue
        lines = body.split('\n')
        while lines and lines[0].startswith('#'):
            h = lines.pop(0)
            blocks.append(dict(t='h2', x=h[3:]) if h.startswith('## ') else dict(t='h1', x=h[2:]))
        if lines: blocks.append(dict(t='p', x=' '.join(lines).strip()))
    return meta, blocks

order = []
def cite(m):
    keys = [k.strip().lstrip('@') for k in m.group(1).split(';')]
    nums = []
    for k in keys:
        if k not in REFS: raise KeyError(k)
        if k not in order: order.append(k)
        nums.append(order.index(k) + 1)
    nums = sorted(set(nums)); parts = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
        parts.append(str(nums[i]) if j == i else (f'{nums[i]}, {nums[j]}' if j == i + 1 else f'{nums[i]}–{nums[j]}')); i = j + 1
    return '(' + ', '.join(parts) + ')'

def sub(x):
    x = re.sub(r'\{\{(\w+)\}\}', lambda m: V[m.group(1)], x)
    return re.sub(r'\[(@[^\]]+)\]', cite, x)

meta, blocks = parse('ms.md')
for b in blocks: b['x'] = sub(b['x'])
_, leg = parse('legends.md')
for b in leg: b['x'] = sub(b['x'])
left = re.findall(r'\{\{\w+\}\}', json.dumps(blocks))
assert not left, left

# Table 1
rows = [['Molecules (K)', 'Exhaustive pairwise screen', 'Combinatorial, size 2', 'Combinatorial, size 1', 'Combinatorial, size 0.5', 'Saving (size 1)']]
for K in (10, 18, 30, 50, 75, 100):
    def f(d):
        v = tp(K, 5, d); return '–' if np.isnan(v) else '%d' % round(v)
    rows.append([str(K), f'{1 + K + comb(K, 2):,}', f(2.0), f(1.0), f(0.5), '%.1f×' % tp(K, 5, 1.0, 'gain')])
table = dict(title='Table 1. Number of mixtures needed to recover 95% of five pairwise interactions',
             note='Interaction size is expressed in replicate standard deviations (Materials and Methods); 3 replicates per mixture; half of the molecules have an effect of their own; diffuse background of small interactions (total variance 0.1). Saving = exhaustive conditions / combinatorial mixtures for size 1. –, 95% recall not reached within 15 K mixtures. The rule N95 ≈ 1.6 K + 6 M + 50/d² − 15 reproduces these values with a median error of 5%.',
             rows=rows)
refs = [REFS[k] for k in order]
json.dump(dict(meta=meta, blocks=blocks, legends=leg, refs=refs, table=table), open('ms.json', 'w'), ensure_ascii=False, indent=1)

def words(t): return len(re.sub(r'[*^~]', '', t).split())
sec = None; wc = {}
for b in blocks:
    if b['t'] == 'h1': sec = b['x']
    elif b['t'] == 'p': wc[sec] = wc.get(sec, 0) + words(b['x'])
print(wc); print('main text (intro+results+discussion):', sum(wc.get(s, 0) for s in ('Introduction', 'Results', 'Discussion')))
print('refs', len(refs)); print(V)
unused = set(REFS) - set(order); print('unused refs', unused)
