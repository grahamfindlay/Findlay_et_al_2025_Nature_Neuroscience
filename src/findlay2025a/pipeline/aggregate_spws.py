from findlay2025a import agg, core, sharp_waves
from findlay2025a.constants import Files


def do_project():
    nb = core.get_project("seahorse")

    spws, _, c_rates = sharp_waves.aggregate_spws()

    c_spws = agg.aggregated_events_wide_to_long(spws)
    c_means = (
        c_spws.groupby(["subject", "experiment", "condition"])
        .mean(numeric_only=True)
        .add_prefix("spw_mean_")
    )
    c_sums = (
        c_spws.groupby(["subject", "experiment", "condition"])
        .sum(numeric_only=True)
        .add_prefix("spw_total_")
    )

    # Prefix for backward compatibility.
    c_rates.rename(
        columns={
            col: f"spw_{col}" for col in c_rates.columns if col not in ["duration"]
        },
        inplace=True,
    )

    contrasts = agg.get_contrasts(
        c_means,
        c_sums,
        c_rates,
        ["spw_count", "spw_rate", "spw_rate_rel2total"],
    )

    spws.to_parquet(nb.get_project_file(Files.CLASSIC_SPWS))
    c_rates.to_parquet(nb.get_project_file(Files.SPW_CONDITION_RATES))
    c_means.to_parquet(nb.get_project_file(Files.SPW_CONDITION_MEANS))
    c_sums.to_parquet(
        nb.get_project_file(Files.SPW_CONDITION_SUMS)
    )  # TODO: Remove, never used.
    contrasts.to_parquet(nb.get_project_file(Files.SPW_CONDITION_CONTRASTS))
