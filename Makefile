.PHONY: test lint typecheck pyright validate synthetic benchmark report fetch tissue all

test:
	pytest

lint:
	ruff format --check .
	ruff check .

typecheck:
	mypy --strict src

pyright:
	pyright

validate:
	alia validate-manifest data/registry

synthetic:
	alia build --synthetic --output data/synthetic_cells.csv

benchmark:
	alia benchmark --input data/synthetic_cells.csv --output results/synthetic

report:
	alia report --results results/synthetic/benchmark.json --output results/synthetic/report.html

fetch:
	alia fetch --registry data/registry --cache-dir $${ALIA_CACHE_DIR:-.cache/alia}

tissue:
	PYTHONPATH=src python scripts/build_tissue_release_assets.py --help

all: validate synthetic benchmark report test
