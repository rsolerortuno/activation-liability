# Critical evaluation — v0.5.0

## 1. Hiring manager, biologics discovery — 9.0/10

**Harshest objection:** the tool now adds epithelial, stromal and endothelial tissue coverage and identifies concrete footprint expansion, but the tissue evidence is RNA-only and therefore cannot establish accessible target density.

**What raises it by one point:** matched target-level surface protein in several paired tissue or stimulation contexts, plus one documented example where the audit changes a modality or target decision.

## 2. Methods reviewer — 9.1/10

**Harshest objection:** the paired design, abstention, marker-leakage controls and count-model sensitivity are defensible, but broad heuristic annotation is weaker than source-author cell labels or an independently validated reference mapping.

**What raises it by one point:** reproduce the tissue result using independent author annotations or a blinded reference mapper and retain the same target ordering and footprint conclusions.

## 3. Skeptical statistician — 8.8/10

**Harshest objection:** there are only five psoriasis pairs, the primary confidence interval uses a normal 1.96 multiplier, and maximum-across-context scoring can favour a single optimistic endpoint.

**What raises it by one point:** Hartung-Knapp or small-sample intervals, hierarchical modelling across contexts, and another untouched paired tissue cohort.

## 4. Senior software engineer — 9.4/10

**Harshest objection:** accession-specific adapters remain procedural and the memory-isolated psoriasis workflow depends on raw GEO layout conventions.

**What raises it by one point:** a checksum-pinned container/lockfile and automated release workflow that rebuilds all public results from a mounted content-addressed cache.

## 5. Reproducibility auditor — 9.3/10

**Harshest objection:** inputs and outputs are hashed and the large archives are reconstructed deterministically, but independent execution on a second physical machine and archived licence evidence are still absent.

**What raises it by one point:** signed second-machine reproduction, container digest, environment lock and archived licence snapshots for every input.

## Weighted overall grade

Weights: biologics discovery 30%, methods 25%, statistics 20%, software 15%, reproducibility 10%.

**Weighted grade: 9.0/10.**

This is a portfolio/research-software grade, not a claim that the method is clinically validated.

## Most likely dismissal reason

A reviewer may still dismiss the project as an RNA-centric tissue-inflammation audit whose most differentiating claim—surface-target liability—has not been confirmed for benchmark targets at protein level.

## What must be true for a 9.5/10

1. Several benchmark controls have matched RNA and surface-protein measurements in paired normal-cell activation or tissue inflammation.
2. An untouched second tissue cohort reproduces target ordering and footprint expansion without marker or threshold changes.
3. A small-sample hierarchical/count model agrees with the paired estimator and gives calibrated intervals.
4. A second machine rebuilds the release from hashes and a locked/containerised environment.

A 10/10 would additionally require prospective or clinically linked evidence connecting inducibility to accessible target density, exposure and observed safety. More transcriptomic cohorts alone cannot justify 10/10.
