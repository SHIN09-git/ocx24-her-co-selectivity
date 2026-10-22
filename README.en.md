# A grouped evaluation of HER premeasurement for CO selectivity prediction in OCx24

[中文](README.md) · [Reproduction](REPRODUCING.md) · [Authorship](AUTHORSHIP.md) · [Licensing](LICENSING.md)

**HER premeasurement did not provide stable incremental predictive value within the fixed OCx24 dataset, inputs and model classes studied here.** This release preserves a negative result and its auditable computational workflow. It does not establish that HER is universally uninformative.

Author: **JIN HANZHE**, Independent Researcher, Chongqing, China. Contact: jinhanzhe09@gmail.com. Original analysis: 2026-10-12. Public packaging: 2026-10-22. **Local v1.0.0 release candidate; this project has not yet been published.** No independent human scientific verification, laboratory replication or peer review is claimed.

The analysis tests whether absolute HER full-cell voltage at 50 mA/cm² improves prediction of CO Faradaic efficiency at 200 mA/cm². M0 uses nominal elemental fractions and synthesis source; M1 adds HER voltage. The primary model is a restricted random forest, with ridge and XRF-based comparisons.

|Primary comparison|Electrodes / preparation IDs / element sets|M0 MAE|M1 MAE|Gain (M0−M1), 95% interval|
|---|---:|---:|---:|---:|
|Development, held-out element sets|404 / 216 / 38|13.76|13.80|−0.04 [−1.79, 1.44]|
|Original held-out evaluation|37 / 14 / 10|11.26|12.33|−1.06 [−3.16, 0.28]|

Errors are preparation-equal MAEs in percentage points. Intervals resample element-set clusters conditional on saved predictions. Values are rounded separately.

## Reproduce

With Python 3.12, the default audit needs only the standard library:

```console
python reproduce.py
```

It verifies release hashes, checks 441 raw pairs and related mapping records, and independently recalculates 23 saved metric tasks without refitting models.

For full data rebuilding and model refitting:

```console
python -m pip install -r requirements.txt
python reproduce.py --full
```

Inputs are bundled. After dependency installation, reproduction is offline. Every run creates a fresh directory under `reproductions/`, or at `--output NEW_DIRECTORY`. Frozen release evidence is never overwritten. See [REPRODUCING.md](REPRODUCING.md) for Windows/macOS/Linux commands and [VALIDATION.md](VALIDATION.md) for the actual tested scope.

## Interpretation and provenance

The 441 electrode pairs represent 230 preparation identifiers and 48 element sets; identifiers do not establish independent syntheses. The split followed a public-data pre-audit. The protocol was locally frozen, not publicly preregistered. The historical filename `preregistered_analysis_plan.json` is preserved as evidence, not a claim of public registration. Replaying already-public outcomes is not a new blind evaluation.

The data cannot establish whether conducting HER causally changes subsequent CO₂ performance. Source holdouts also change elemental coverage. The ethylene-rate appendix concerns only a pinned historical code version and molar rates, not a sixfold error in ethylene FE or an explanation of the primary CO result.

Data: [Abed et al., OCx24](https://arxiv.org/abs/2411.11783), pinned FAIR-Chem commit `0a625a8924cf0a87dffbcd2f0ac58d82ca43bfa4`. [Official dataset license](https://facebookresearch.github.io/fairchem/ocx24/): CC BY 4.0. Original project code: MIT. Original documentation, figures and derived tables: CC BY 4.0. Bundled upstream code retains Meta's MIT notice. Cite both this analysis and OCx24; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and [CITATION.cff](CITATION.cff).

Analysis, source checks, drafting and packaging used OpenAI Codex. AI is not an author; see [AUTHORSHIP.md](AUTHORSHIP.md).

