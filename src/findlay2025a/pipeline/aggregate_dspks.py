from findlay2025a import agg, core, dentate_spikes
from findlay2025a.constants import Files


def do_project():
    nb = core.get_project("seahorse")

    dspks, _, _, c_rates = dentate_spikes.aggregate_dspks()
    c_dspks = agg.aggregated_events_wide_to_long(dspks)
    c_means = (
        c_dspks.groupby(["subject", "experiment", "condition"])
        .mean(numeric_only=True)
        .add_prefix("dspk_mean_")
    )
    c_sums = (
        c_dspks.groupby(["subject", "experiment", "condition"])
        .sum(numeric_only=True)
        .add_prefix("dspk_total_")
    )

    # Prefix for backward compatibility.
    c_rates.rename(
        columns={
            col: f"dspk_{col}" for col in c_rates.columns if col not in ["duration"]
        },
        inplace=True,
    )

    contrasts = agg.get_contrasts(
        c_means,
        c_sums,
        c_rates,
        ["dspk_count", "dspk_rate", "dspk_rate_rel2total"],
    )

    dspks.to_parquet(nb.get_project_file(Files.DSPKS))
    c_rates.to_parquet(nb.get_project_file(Files.DSPK_CONDITION_RATES))
    c_means.to_parquet(nb.get_project_file(Files.DSPK_CONDITION_MEANS))
    contrasts.to_parquet(nb.get_project_file(Files.DSPK_CONDITION_CONTRASTS))
