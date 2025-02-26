from findlay2025a import agg, core
from findlay2025a.constants import Files


def do_project():
    nb = core.get_project("seahorse")

    pwr = agg.aggregate_cortical_bandpowers()
    c_pwr = agg.aggregated_events_wide_to_long(pwr)
    c_means = (
        c_pwr.groupby(["subject", "experiment", "condition"])
        .mean(numeric_only=True)
        .add_prefix("cx_mean_")
    )
    c_sums = (
        c_pwr.groupby(["subject", "experiment", "condition"])
        .sum(numeric_only=True)
        .add_prefix("cx_total_")
    )
    contrasts = agg.get_contrasts(c_means)

    c_means.to_parquet(nb.get_project_file(Files.CX_BANDPOWER_MEANS))
    c_sums.to_parquet(
        nb.get_project_file(Files.CX_BANDPOWER_SUM)
    )  # TODO: Remove, never used.
    contrasts.to_parquet(nb.get_project_file(Files.CX_BANDPOWER_CONTRASTS))
