from findlay2025a import agg, core
from findlay2025a.constants import Files


def do_project():
    nb = core.get_project("seahorse")

    emg = agg.aggregate_emg()
    c_pwr = agg.aggregated_events_wide_to_long(emg)
    c_means = (
        c_pwr.groupby(["subject", "experiment", "condition"])
        .mean(numeric_only=True)
        .add_prefix("mean_")
    )
    c_sums = (
        c_pwr.groupby(["subject", "experiment", "condition"])
        .sum(numeric_only=True)
        .add_prefix("total_")
    )

    c_means.to_parquet(nb.get_project_file(Files.EMG_CONDITION_MEANS))
    c_sums.to_parquet(nb.get_project_file(Files.EMG_CONDITION_SUMS))
