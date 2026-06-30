.PHONY: setup data run analyze crosscheck all clean

CONDA_ENV ?= epi-clock-replication
PYTHON ?= python

setup:
	conda env create -f environment.yml || conda env update -f environment.yml --prune

data:
	$(PYTHON) src/download_data.py

analyze:
	$(PYTHON) src/analyze.py

crosscheck:
	$(PYTHON) src/crosscheck.py

run: analyze crosscheck

all: data run

clean:
	$(PYTHON) -c "import shutil, pathlib; results=pathlib.Path('results'); shutil.rmtree(results, ignore_errors=True); results.mkdir(exist_ok=True); (results/'.gitkeep').touch(); [p.unlink(missing_ok=True) for p in [pathlib.Path('data/raw/betas_subset.parquet'), pathlib.Path('data/raw/samples_metadata.parquet')]]"
