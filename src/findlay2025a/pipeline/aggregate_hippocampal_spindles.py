from findlay2025a import agg, constants, core, spindles


def do_project():
    nb = core.get_project("seahorse")

    sps, _, _, c_rates = spindles.aggregate_spindles(region="hippocampal")

    c_sps = agg.aggregated_events_wide_to_long(sps)
    c_means = (
        c_sps.groupby(["subject", "experiment", "condition"])
        .mean(numeric_only=True)
        .add_prefix("hcsp_mean_")
    )
    c_sums = (
        c_sps.groupby(["subject", "experiment", "condition"])
        .sum(numeric_only=True)
        .add_prefix("hcsp_total_")
    )

    # Prefix for backward compatibility.
    c_rates.rename(
        columns={
            col: f"hcsp_{col}" for col in c_rates.columns if col not in ["duration"]
        },
        inplace=True,
    )

    contrasts = agg.get_contrasts(
        c_means,
        c_sums,
        c_rates,
        ["hcsp_count", "hcsp_rate", "hcsp_rate_rel2total"],
    )

    sps.to_parquet(nb.get_project_file(constants.Files.HIPPOCAMPAL_SPINDLES))
    c_rates.to_parquet(
        nb.get_project_file(constants.Files.HIPPOCAMPAL_SPINDLE_CONDITION_RATES)
    )
    c_means.to_parquet(
        nb.get_project_file(constants.Files.HIPPOCAMPAL_SPINDLE_CONDITION_MEANS)
    )
    contrasts.to_parquet(
        nb.get_project_file(constants.Files.HIPPOCAMPAL_SPINDLE_CONDITION_CONTRASTS)
    )
