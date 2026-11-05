#!/usr/bin/env python3
"""Rebuild the annual service and biochemical-threshold figures."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "aggregate"
OUTPUT = ROOT / "figures"
OUTPUT.mkdir(exist_ok=True)

GENOTYPES = [
    "No mutations",
    "H63D carrier",
    "H63D homozygote",
    "C282Y carrier",
    "Compound heterozygote",
    "C282Y Homozygous",
]
COLOURS = {
    "No mutations": "#D0D0D0",
    "H63D carrier": "#808080",
    "H63D homozygote": "#2C7FB8",
    "C282Y carrier": "#404040",
    "Compound heterozygote": "#7B3294",
    "C282Y Homozygous": "#D7191C",
}


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUTPUT / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def annual_service_figure() -> None:
    distribution = pd.read_csv(DATA / "annual_genotype_distribution.csv")
    positivity = pd.read_csv(DATA / "annual_genotype_positivity.csv")

    years = positivity["genetic_test_year"].to_numpy()
    percentages = distribution.pivot(
        index="genetic_test_year", columns="genotype", values="genotype_percent"
    ).reindex(index=years, columns=GENOTYPES)

    fig, axis = plt.subplots(figsize=(15, 7.2), facecolor="white")
    bottom = np.zeros(len(years))
    for genotype in GENOTYPES:
        values = percentages[genotype].to_numpy()
        axis.bar(
            years,
            values,
            bottom=bottom,
            color=COLOURS[genotype],
            edgecolor="white",
            linewidth=0.5,
            width=0.82,
        )
        bottom += values

    line_specs = [
        ("positivity_percent", "Genotype positivity", "#000000"),
        ("c282y_homozygous_percent", "C282Y/C282Y", "#D7191C"),
        ("compound_heterozygote_percent", "C282Y/H63D", "#7B3294"),
        ("h63d_homozygote_percent", "H63D/H63D", "#2C7FB8"),
    ]
    for column, _, colour in line_specs:
        axis.plot(
            years,
            positivity[column],
            color=colour,
            marker="o",
            linewidth=2.4,
            markersize=4.5,
        )

    count_axis = axis.twinx()
    count_axis.plot(
        years,
        positivity["annual_genotyped_patients"],
        color="#F18701",
        marker="^",
        linewidth=2.4,
        markersize=5,
    )

    axis.set(
        title="Annual distribution of HFE genotypes and positivity, 2004-2024",
        xlabel="Genetic test year",
        ylabel="Percentage of genotyped patients",
        xlim=(years.min() - 0.7, years.max() + 0.7),
        ylim=(0, 102),
    )
    count_axis.set_ylabel("First-tested patients per year")
    axis.set_xticks(years)
    axis.tick_params(axis="x", rotation=45)
    axis.grid(axis="y", color="#E5E7E9", linewidth=0.7)
    axis.set_axisbelow(True)

    handles = [
        Line2D([0], [0], color=colour, marker="o", linewidth=2, label=label)
        for _, label, colour in line_specs
    ]
    handles.append(
        Line2D([0], [0], color="#F18701", marker="^", linewidth=2, label="First-tested patients")
    )
    handles.extend(
        Patch(facecolor=COLOURS[genotype], label=genotype) for genotype in reversed(GENOTYPES)
    )
    axis.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=4,
        frameon=False,
        fontsize=8.5,
    )
    fig.subplots_adjust(bottom=0.28)
    save(fig, "annual_genotype_distribution_and_positivity")


def threshold_figures() -> None:
    composition = pd.read_csv(DATA / "easl_genotype_composition.csv")
    status_order = list(dict.fromkeys(composition["threshold_status"]))

    fig, axis = plt.subplots(figsize=(8.2, 4.8), facecolor="white")
    left = np.zeros(len(status_order))
    for genotype in GENOTYPES:
        rows = (
            composition.loc[composition["genotype"].eq(genotype)]
            .set_index("threshold_status")
            .reindex(status_order)
        )
        values = rows["genotype_percent"].to_numpy()
        axis.barh(
            status_order,
            values,
            left=left,
            color=COLOURS[genotype],
            edgecolor="white",
            label=genotype,
        )
        left += values
    axis.set(xlabel="Percentage of patients", xlim=(0, 100))
    axis.set_title("Genotype composition by biochemical-threshold status")
    axis.grid(axis="x", color="#E5E7E9", linewidth=0.7)
    axis.set_axisbelow(True)
    axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False)
    fig.subplots_adjust(bottom=0.27, left=0.27)
    save(fig, "easl_genotype_composition")

    annual = pd.read_csv(DATA / "easl_annual_summary.csv")
    fig, axis = plt.subplots(figsize=(9.2, 4.8), facecolor="white")
    axis.plot(
        annual["genetic_test_year"],
        annual["threshold_concordance_percent"],
        color="#2166AC",
        marker="o",
        label="Meeting both thresholds",
    )
    axis.plot(
        annual["genetic_test_year"],
        annual["positivity_threshold_positive_percent"],
        color="#B2182B",
        marker="o",
        label="Genetic positivity: both thresholds",
    )
    axis.plot(
        annual["genetic_test_year"],
        annual["positivity_threshold_negative_percent"],
        color="#666666",
        marker="o",
        label="Genetic positivity: not both",
    )
    axis.set(
        title="Biochemical-threshold attainment and genetic yield, 2004-2024",
        xlabel="Genetic test year",
        ylabel="Patients (%)",
        ylim=(0, 100),
    )
    axis.set_xticks(annual["genetic_test_year"])
    axis.tick_params(axis="x", rotation=45)
    axis.grid(color="#E5E7E9", linewidth=0.7)
    axis.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=3)
    fig.subplots_adjust(bottom=0.30)
    save(fig, "easl_threshold_attainment_and_yield")


if __name__ == "__main__":
    annual_service_figure()
    threshold_figures()
    print(f"Figures written to {OUTPUT}")
