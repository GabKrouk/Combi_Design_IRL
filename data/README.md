# Raw data (not versioned)

Place the public files below in this folder, then run `python prepare_data.py` from `analysis/`.

```
data/
  ristova/GSE71737_series_matrix.txt.gz
      GEO series matrix of GSE71737 (Ristova et al. 2016, Sci. Signal. 9:rs13)
      https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71737
  kellermeier/Kellermeier2014_Supplemental_DS1.xlsx   (root traits, 2^5 N/P/K/S/light)
  kellermeier/Kellermeier2014_Supplemental_DS4.xlsx   (transcriptome, 2^3 N/P/K)
      Supplemental Data Sets of Kellermeier et al. 2014, Plant Cell 26:1480-1496
  tekin/Tekin2018_Dataset_MOESM2.pdf
      Supplementary dataset of Tekin et al. 2018, NPJ Syst. Biol. Appl. 4:31
      (growth for every combination within each five-antibiotic panel; parsed with `pdftotext -raw`)
```

`pdftotext` (poppler-utils) is needed to read the Tekin dataset.
