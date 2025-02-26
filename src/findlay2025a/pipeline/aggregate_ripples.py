from findlay2025a import agg, core, ripples
from findlay2025a.constants import Files


def do_project():
    nb = core.get_project("seahorse")

    rips, _, c_rates = ripples.aggregate_ripples()
    c_rips = agg.aggregated_events_wide_to_long(rips)
    c_means = (
        c_rips.groupby(["subject", "experiment", "condition"])
        .mean(numeric_only=True)
        .add_prefix("ripple_mean_")
    )
    c_sums = (
        c_rips.groupby(["subject", "experiment", "condition"])
        .sum(numeric_only=True)
        .add_prefix("ripple_total_")
    )

    # Prefix for backward compatibility.
    c_rates.rename(
        columns={
            col: f"ripple_{col}" for col in c_rates.columns if col not in ["duration"]
        },
        inplace=True,
    )

    contrasts = agg.get_contrasts(
        c_means,
        c_sums,
        c_rates,
        ["ripple_count", "ripple_rate", "ripple_rate_rel2total"],
    )

    rips.to_parquet(nb.get_project_file(Files.RIPPLES))
    c_rates.to_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_RATES))
    c_means.to_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_MEANS))
    c_sums.to_parquet(
        nb.get_project_file(Files.RIPPLE_CONDITION_SUMS)
    )  # TODO: Remove, never used.
    contrasts.to_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_CONTRASTS))
