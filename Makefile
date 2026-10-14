PY = python3
A = cd analysis &&

.PHONY: all data bench sims figures manuscript

all: data bench sims figures manuscript

data:
	$(A) $(PY) prepare_data.py && $(PY) truth.py

bench:
	$(A) $(PY) run_bench.py ristova res_ristova.pkl
	$(A) $(PY) run_bench.py tekin res_tekin.pkl
	$(A) $(PY) run_bench.py kellermeier res_kellermeier.pkl
	$(A) $(PY) figures.py && $(PY) by_order.py && $(PY) rarity.py && $(PY) eb.py && $(PY) recal_check.py

sims:
	$(A) $(PY) simK.py trip fig/simK_trip.csv 2
	$(A) $(PY) simK2.py pairs fig/simK_pairs.csv
	$(A) $(PY) simK2.py bg fig/simK_bg.csv
	$(A) $(PY) simK_ex.py

figures:
	$(A) $(PY) paper_figs.py

manuscript:
	cd manuscript && $(PY) fill.py && node build_docx.js && node build_pptx.js
