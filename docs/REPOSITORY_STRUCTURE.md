# Repository structure and reproducibility contract

`activation-liability` is distributed as an installable Python application. The GitHub repository
contains everything needed to inspect, test and run the software, except the multi-gigabyte raw
public datasets.

```text
activation-liability/
├── .github/workflows/ci.yml       # real CI: lint, typing, CLI smoke tests and pytest
├── config/                         # frozen benchmark/evidence configuration
├── data/
│   ├── candidates/                 # source catalog and next-phase expectations
│   ├── controls/                   # predeclared positive/negative controls
│   └── registry/                   # dataset manifests and validation schemas
├── docs/                           # methods, assumptions, decisions and limitations
│   └── assets/                     # versioned SVG figures used by the README
├── notebooks/                      # resumable public-data download notebook
├── results/
│   ├── synthetic/                  # small reproducible demonstration outputs
│   └── real/                       # compact executed benchmark artefacts
├── scripts/                        # data bootstrap and release-building utilities
├── src/activation_liability/       # importable package and `alia` CLI
├── tests/                          # 52 automated tests
├── pyproject.toml                  # dependencies and console-script entry point
├── LICENSE                         # Apache-2.0
└── README.md
```

## What a user can run immediately

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -e '.[dev]'
alia --help
pytest -q
```

A synthetic end-to-end run requires no external data:

```bash
alia build --synthetic --output /tmp/alia_synthetic.csv
alia audit --input /tmp/alia_synthetic.csv --output results/local_audit
alia benchmark --input /tmp/alia_synthetic.csv --output results/local_benchmark
alia report --results results/local_benchmark/benchmark.json --output results/local_report.html
```

## What is not committed

The repository does not include raw GEO TAR files, FASTQ/BAM files, large MatrixMarket inputs,
virtual environments or credentials. This is intentional. Official download URLs, checksums,
Drive locations and ingestion assumptions are recorded in `docs/DATA.md`,
`docs/DRIVE_DATA_INDEX.md`, the dataset registry and the committed result manifests.

## Reproducibility levels

1. **Software reproduction:** install the package and run all tests and synthetic examples.
2. **Result inspection:** inspect committed compact tables, JSON claims, reports and checksums.
3. **Real-data reconstruction:** download the named public files and run `alia real-validate` and
   `alia tissue-validate` according to the runbook.
4. **External reproduction:** run the same release on a clean second machine and compare output
   manifests.
