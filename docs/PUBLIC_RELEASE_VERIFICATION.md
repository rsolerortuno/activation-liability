# Public release verification — 2026-07-31

The public package was validated after adding documentation, graphics and the GSE190564 download
contract.

```text
ruff format --check .                         PASS
ruff check .                                  PASS
mypy --strict src                             PASS — 21 modules
Pyright 1.1.411                               PASS — 0 errors, 0 warnings
pytest                                        PASS — 52 tests
coverage (configured non-I/O scope)           PASS — 85.38%
README quantitative truthfulness tests         PASS
GSE190564 required-modality validator tests    PASS
```

No GSE190564 data or result is included. The notebook and validator only establish download
integrity and required GEX/ADT presence.
