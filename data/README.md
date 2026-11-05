# Aggregate data

These CSV files contain non-identifying results used to reconstruct selected poster figures.

| File | Content |
|---|---|
| `annual_genotype_distribution.csv` | Annual patient counts by reportable genotype, 2004-2024 |
| `annual_genotype_positivity.csv` | Annual testing volume and positive-genotype counts |
| `easl_genotype_composition.csv` | Genotype counts by biochemical-threshold status |
| `easl_annual_summary.csv` | Annual threshold attainment and genetic yield |
| `genotype_roc_metrics.csv` | Out-of-fold AUC summaries for the four genotype outcomes |
| `genotype_roc_coordinates.csv` | ROC curves interpolated to a fixed false-positive-rate grid |
| `biochemical_model_metrics.csv` | Complete-profile model performance |
| `biochemical_roc_coordinates.csv` | Complete-profile ROC curves on a fixed grid |
| `biochemical_feature_importance.csv` | Aggregate permutation-importance estimates |
| `biochemical_effects.csv` | Standardised biochemical group differences |
| `referral_feature_associations.csv` | Adjusted odds ratios and confidence intervals |

Genetic positivity comprises C282Y/C282Y, C282Y/H63D and H63D/H63D. The underlying cohort construction, deduplication, linkage and model fitting were performed on secure patient-level data and cannot be reconstructed from these files.
