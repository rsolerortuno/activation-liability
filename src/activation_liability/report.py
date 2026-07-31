"""Dependency-free HTML result renderer."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any


def render_report(payload: dict[str, Any], output: Path) -> None:
    """Render a compact HTML report that remains honest about result status."""

    title = "activation-liability report"
    status = html.escape(str(payload.get("status", "UNKNOWN")))
    if "target_summary" in payload:
        rows = payload["target_summary"]
    else:
        rows = payload.get("target_scores", [])
    table_rows = []
    for row in rows[:25]:
        table_rows.append(
            "<tr>"
            f"<td>{html.escape(str(row.get('target', '')))}</td>"
            f"<td>{html.escape(str(round(float(row.get('score', 0.0)), 3)))}</td>"
            f"<td>{html.escape(str(row.get('evidence_class', row.get('tier', ''))))}</td>"
            f"<td>{html.escape(str(row.get('protein_concordance', '')))}</td>"
            "</tr>"
        )
    raw = html.escape(json.dumps(payload, indent=2, sort_keys=True))
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body {{font-family:system-ui;max-width:1100px;margin:40px auto;padding:0 20px;line-height:1.5}}
table {{border-collapse:collapse;width:100%}}
th,td {{border:1px solid #ddd;padding:8px;text-align:left}}
.warning {{padding:12px;border:1px solid #b66;background:#fff7f2}}
pre {{overflow:auto;background:#f6f6f6;padding:16px}}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="warning">
Status: <strong>{status}</strong>. This is a liability detector, not a toxicity predictor.
</p>
<table>
<thead><tr><th>Target</th><th>Score</th><th>Evidence/tier</th><th>Protein</th></tr></thead>
<tbody>{"".join(table_rows)}</tbody>
</table>
<h2>Machine-readable payload</h2>
<pre>{raw}</pre>
</body>
</html>"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")
