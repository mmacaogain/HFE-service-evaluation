# HFE service evaluation

Preliminary materials accompanying the [ACBI 2026 poster](poster/ACBI_2026_HFE_service_evaluation.pdf) and [conference abstract](abstract/Mac%20Aog%C3%A1in%20et%20al.%202026%20ACBI%20abstract.docx) evaluating HFE genetic testing at the St James's Hospital Biochemistry Department for the 48th Annual Conference.

This repository contains the poster, disclosure-checked aggregate results and code that rebuilds selected figures. It does not contain patient-level data, the production linkage pipeline or a validated clinical prediction tool.

## Contents

- `poster/`: [conference poster](poster/ACBI_2026_HFE_service_evaluation.pdf)
- `abstract/`: [conference abstract](abstract/Mac%20Aog%C3%A1in%20et%20al.%202026%20ACBI%20abstract.docx)
- `data/aggregate/`: non-identifiable aggregate inputs used by the plotting code
- `analysis/`: figure-reproduction scripts
- `figures/`: regenerated outputs

## Reproduce

```bash
python analysis/01_service_and_thresholds.py
python analysis/02_models_and_referrals.py
```

Install the Python packages listed in `requirements.txt` first. The ferritin-TSAT patient scatterplot is retained in the poster but is not regenerated because it requires patient-level measurements.

## Status and data availability

These are preliminary poster-reproduction materials primarily intended for conference delegates who wish to better understand or review the analysis underpinning the poster presentation. Clinical source data and pseudonymised patient-level records cannot be shared at this time. The released aggregate files support checking the plotted values and figure construction, but not refitting the underlying patient-level models.

OpenAI Codex assisted with code development and visualisation. Identifiable data remained within the secure institutional environment.
