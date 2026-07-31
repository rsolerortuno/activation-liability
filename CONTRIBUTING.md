# Contributing

Contributions are welcome when they improve reproducibility, statistical validity, data provenance
or scientific clarity.

1. Open an issue describing the proposed change.
2. Keep raw and controlled-access data out of Git.
3. Add or update tests for software changes.
4. Run `ruff format --check .`, `ruff check .`, `mypy --strict`, `pyright` and `pytest`.
5. Document any change to evidence rules, thresholds or permitted claims.
6. Never retune a frozen benchmark after inspecting a new holdout cohort without labelling the
   result exploratory.

Dataset proposals should include accession, modality, donor/patient pairing, reuse terms and the
specific evidence gap they could close.
