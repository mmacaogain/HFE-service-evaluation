#!/usr/bin/env python3
"""Rebuild the ROC, biochemical-profile and referral-feature figures."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "aggregate"
OUTPUT = ROOT / "figures"
OUTPUT.mkdir(exist_ok=True)

ROC_OUTCOMES = [
    ("all_positive", "A  All positive"),
    ("c282y_homozygous_with_carriers", "B  C282Y/C282Y"),
    ("compound_heterozygous", "C  C282Y/H63D"),
    ("h63d_homozygous", "D  H63D/H63D"),
]
ROC_MODELS = ["All five combined", "TSAT", "Ferritin", "Family history", "Age", "Sex"]
ROC_COLOURS = {
    "All five combined": "#111111",
    "TSAT": "#D7191C",
    "Ferritin": "#7B3294",
    "Family history": "#009999",
    "Age": "#2C7BB6",
    "Sex": "#6B6B6B",
}
ROC_STYLES = {"Family history": "--", "Sex": "--"}


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUTPUT / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def genotype_roc_figure() -> None:
    curves = pd.read_csv(DATA / "genotype_roc_coordinates.csv")
    metrics = pd.read_csv(DATA / "genotype_roc_metrics.csv")

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.6), sharex=True, sharey=True)
    for axis, (outcome, title) in zip(axes.flat, ROC_OUTCOMES):
        outcome_metrics = metrics.loc[metrics["outcome"].eq(outcome)].set_index("model")
        for model in ROC_MODELS:
            points = curves.loc[
                curves["outcome"].eq(outcome) & curves["model"].eq(model)
            ]
            axis.plot(
                points["false_positive_rate"],
                points["true_positive_rate"],
                color=ROC_COLOURS[model],
                linestyle=ROC_STYLES.get(model, "-"),
                linewidth=2.8 if model == "All five combined" else 1.8,
            )
        combined = outcome_metrics.loc["All five combined"]
        axis.plot([0, 1], [0, 1], color="#B9B9B9", linestyle=":")
        axis.set_title(title, loc="left", weight="bold", fontsize=14, pad=18)
        axis.text(
            0,
            1.025,
            f"n={int(combined.patients):,}; cases={int(combined.positive_patients):,}",
            transform=axis.transAxes,
            fontsize=9,
            color="#55616D",
        )
        axis.text(
            0.96,
            0.05,
            f"Combined AUC {combined.oof_auc:.3f}\n95% CI {combined.ci_95_low:.3f}-{combined.ci_95_high:.3f}",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.7,
            bbox={"facecolor": "white", "edgecolor": "#D9D9D9", "alpha": 0.9},
        )
        axis.set(xlim=(0, 1), ylim=(0, 1.01))
        axis.grid(color="#E5E7E9", linewidth=0.7)
        axis.set_aspect("equal", adjustable="box")

    axes[0, 0].set_ylabel("True-positive rate")
    axes[1, 0].set_ylabel("True-positive rate")
    axes[1, 0].set_xlabel("False-positive rate")
    axes[1, 1].set_xlabel("False-positive rate")
    handles = [
        Line2D(
            [0],
            [0],
            color=ROC_COLOURS[model],
            linestyle=ROC_STYLES.get(model, "-"),
            linewidth=2.5,
            label="Combined model" if model == "All five combined" else model,
        )
        for model in ROC_MODELS
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
    fig.subplots_adjust(bottom=0.16, hspace=0.30, wspace=0.18)
    save(fig, "genotype_roc_panels")


def biochemical_profile_figure() -> None:
    curves = pd.read_csv(DATA / "biochemical_roc_coordinates.csv")
    metrics = pd.read_csv(DATA / "biochemical_model_metrics.csv").set_index("model")
    importance = pd.read_csv(DATA / "biochemical_feature_importance.csv")
    effects = pd.read_csv(DATA / "biochemical_effects.csv")

    model_order = ["reference", "elastic_net", "random_forest"]
    model_colours = {"reference": "#666666", "elastic_net": "#2166AC", "random_forest": "#D95F02"}
    model_labels = {
        "reference": "Ferritin + TSAT",
        "elastic_net": "9-analyte elastic net",
        "random_forest": "9-analyte random forest",
    }

    order = importance.sort_values("mean_auc_decrease", ascending=False)["analyte"].tolist()
    importance = importance.set_index("analyte").loc[order]
    effects = effects.set_index("analyte").loc[order]
    positions = np.arange(len(order))

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.8), gridspec_kw={"width_ratios": [1, 1.25, 1.15]})
    roc_axis, importance_axis, effect_axis = axes

    for model in model_order:
        points = curves.loc[curves["model"].eq(model)]
        row = metrics.loc[model]
        roc_axis.plot(
            points["false_positive_rate"],
            points["true_positive_rate"],
            color=model_colours[model],
            linewidth=2.3,
            linestyle="--" if model == "reference" else "-",
            label=f"{model_labels[model]}: {row.roc_auc:.3f}",
        )
    roc_axis.plot([0, 1], [0, 1], color="#B9B9B9", linestyle=":")
    roc_axis.set(
        title="A  Out-of-fold ROC",
        xlabel="False-positive rate",
        ylabel="True-positive rate",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    roc_axis.grid(color="#E5E7E9", linewidth=0.7)
    roc_axis.legend(frameon=False, fontsize=8)

    values = importance["mean_auc_decrease"].to_numpy()
    low = importance["empirical_95_low"].to_numpy()
    high = importance["empirical_95_high"].to_numpy()
    importance_axis.errorbar(
        values,
        positions,
        xerr=np.vstack([values - low, high - values]),
        fmt="o",
        color="#D95F02",
        capsize=3,
    )
    importance_axis.axvline(0, color="#888888", linestyle="--")
    importance_axis.set_yticks(positions, importance["analyte_label"])
    importance_axis.invert_yaxis()
    importance_axis.set(title="B  Predictive contribution", xlabel="Held-out AUC decrease")
    importance_axis.grid(axis="x", color="#E5E7E9", linewidth=0.7)

    effect_values = effects["hedges_g"].to_numpy()
    effect_low = effects["ci_95_low"].to_numpy()
    effect_high = effects["ci_95_high"].to_numpy()
    for position, value, low_value, high_value in zip(
        positions, effect_values, effect_low, effect_high
    ):
        colour = "#B2182B" if value >= 0 else "#2166AC"
        effect_axis.errorbar(
            value,
            position,
            xerr=np.array([[value - low_value], [high_value - value]]),
            fmt="o",
            color=colour,
            capsize=3,
        )
    effect_axis.axvline(0, color="#888888", linestyle="--")
    effect_axis.set_yticks(positions, importance["analyte_label"])
    effect_axis.invert_yaxis()
    effect_axis.set(title="C  Direction of difference", xlabel="Hedges' g")
    effect_axis.grid(axis="x", color="#E5E7E9", linewidth=0.7)

    fig.tight_layout()
    save(fig, "c282y_biochemical_profile")


def referral_forest_figure() -> None:
    data = pd.read_csv(DATA / "referral_feature_associations.csv")
    positions = np.arange(len(data))[::-1]

    fig, axis = plt.subplots(figsize=(7.4, 4.7), facecolor="white")
    for position, row in zip(positions, data.itertuples(index=False)):
        colour = "#7B3294" if bool(row.significant) else "#666666"
        axis.errorbar(
            row.adjusted_or,
            position,
            xerr=np.array([[row.adjusted_or - row.ci_low], [row.ci_high - row.adjusted_or]]),
            fmt="o",
            color=colour,
            capsize=3,
        )
        axis.text(row.ci_high * 1.08, position, row.estimate_label, va="center", fontsize=8.5)

    axis.axvline(1, color="#888888", linestyle="--")
    axis.set_xscale("log")
    axis.set_xlim(0.18, 8)
    axis.set_xticks([0.25, 0.5, 1, 2, 4], ["0.25", "0.50", "1.00", "2.00", "4.00"])
    axis.set_yticks(positions, data["feature"])
    axis.set_xlabel("Adjusted odds ratio for C282Y/C282Y (log scale)")
    axis.grid(axis="x", color="#E5E7E9", linewidth=0.7)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)
    fig.tight_layout()
    save(fig, "referral_feature_forest_plot")


if __name__ == "__main__":
    genotype_roc_figure()
    biochemical_profile_figure()
    referral_forest_figure()
    print(f"Figures written to {OUTPUT}")
