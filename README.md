# Combi_Design_IRL

**Finding interactions among many molecules with a small fraction of exhaustive experiments**
Gabriel Krouk & Dennis E. Shasha

This repository contains everything behind the manuscript: the combinatorial designs, the simulations
(up to 100 molecules), the benchmark on three real complete factorial experiments ("in real life", hence *IRL*),
the calibration of the planning model, the figures, the manuscript sources, and a small web calculator.

## The idea in one paragraph

Testing every combination of *K* molecules needs 2^K mixtures; testing every pair still needs
1 + K + K(K−1)/2 conditions (5,051 for 100 molecules). In a **combinatorial design**, each molecule is present
in half of the mixtures and every pair of molecules appears in all four presence/absence states, so every
mixture informs every molecule and every pair (Shasha *et al.*, *Plant Physiol.* 2001). Analyzed with a
parsimonious model (all single-molecule effects + forward selection of interactions), such designs find
interactions with a number of mixtures that grows about linearly with *K*:

```
N95 ≈ 1.6 K + 6 M + 50 / d² − 15
```

*N95* = mixtures (3 replicates each) needed to recover 95% of *M* pairwise interactions of size *d*
(in replicate standard deviations). See `manuscript/` for details and caveats (weak interactions,
dense backgrounds of small interactions, three-way interactions).

## Repository layout

```
analysis/            Python code (run from this folder)
  lib.py               model terms, truth fit, forward-selection analysis, scoring
  designs.py           covering arrays, combinatorial / random / one-at-a-time sequences
  datasets.py          loaders for the three benchmark datasets
  prepare_data.py      raw downloads (data/) -> analysis/prep/*.pkl
  truth.py             reference answers from the complete factorials
  run_bench.py         benchmark of sub-designs on real data (res_*.pkl)
  rarity.py            weight of each interaction order (Fig. 4)
  eb.py                empirical-Bayes effect-size mixtures (calibration)
  recal_check.py       calibrated vs naive simulation vs reality (Fig. 5A-B)
  simK.py / simK2.py   planning simulations, K = 10..100 (pairs) and 10..30 (triplets)
  simK_ex.py           exhaustive pairwise screen with the same analysis
  sim18.py             18-molecule planning simulation
  paper_figs.py        all manuscript figures (paper/*.png, *.pdf) and threshold tables
  figures.py, by_order.py   benchmark summary tables
  fig/                 result tables (CSV) produced by the scripts
  paper/               final figures (PNG 600 dpi + PDF) and tables
  res_*.pkl            benchmark outputs (so figures can be rebuilt without rerunning)
data/                raw public data (not versioned, see data/README.md)
manuscript/          ms.md + legends.md (sources), fill.py, build_docx.js, build_pptx.js,
                     Krouk_Shasha_manuscript.docx, Krouk_Shasha_figures_A4.pptx
web/                 interactive calculator (open web/index.html)
```

## Reproducing

```bash
pip install -r requirements.txt          # numpy, scipy, pandas, scikit-learn, matplotlib, openpyxl
npm install                              # docx, pptxgenjs (manuscript build only)

# 1. raw data -> see data/README.md, then
cd analysis
python prepare_data.py                   # builds prep/*.pkl
python truth.py                          # reference answers
python run_bench.py ristova res_ristova.pkl
python run_bench.py tekin res_tekin.pkl
python run_bench.py kellermeier res_kellermeier.pkl
python figures.py && python by_order.py && python rarity.py && python eb.py && python recal_check.py
python simK.py trip fig/simK_trip.csv 2
python simK2.py pairs fig/simK_pairs.csv && python simK2.py bg fig/simK_bg.csv
python simK_ex.py
python paper_figs.py                     # all figures -> paper/

# 2. manuscript
cd ../manuscript
python fill.py && node build_docx.js && node build_pptx.js
```

Or simply `make all` (see `Makefile`). The simulations take about one hour on two cores; the
benchmark about one hour. Steps after `prepare_data.py` only need the files already in the repository
if you skip the expensive runs (`make figures manuscript`).

## Data sources

| Dataset | System | Factors | Outputs | Source |
|---|---|---|---|---|
| Ristova *et al.* 2016, *Sci. Signal.* | *Arabidopsis* root transcriptome | nitrate, ammonium, auxin, cytokinin, ABA (2^5) | 15,348 genes | GEO GSE71737 |
| Kellermeier *et al.* 2014, *Plant Cell* | *Arabidopsis* root architecture | N, P, K, S, light (2^5) | 15 traits | Supplemental Data Sets 1 and 4 |
| Tekin *et al.* 2018, *NPJ Syst. Biol. Appl.* | *E. coli* growth | 5 of 8 antibiotics (2^5 per panel) | 10,891 panels | Supplementary data |

## Citation

If you use the designs or the calculator, please cite Shasha *et al.* (2001) *Plant Physiol.* 127:1590–1594
(doi:10.1104/pp.010683) and the manuscript (see `CITATION.cff`).

## License

Code: MIT (see `LICENSE`). Data belong to their original authors.
