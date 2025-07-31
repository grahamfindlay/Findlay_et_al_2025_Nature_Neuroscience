import numpy as np
import scipy.stats

from ecephys import wne
from findlay2025a import core, hypnograms
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    hg = hypnograms.load_consolidated_hypnogram(
        experiment, sglx_subject.name, simplify=True, clean=True
    )
    rips = core.read_ripples(sglx_subject.name, experiment, kind="raw")
    rips["state"] = hg.get_states(rips["pk_time"])
    rips = rips[rips["duration"] <= 0.350].reset_index(drop=True)
    rips = rips[rips["pk_amp"] <= 800].reset_index(drop=True)
    rips["zlog_amp"] = scipy.stats.zscore(np.log(rips["pk_amp"]))
    rips["zlog_duration"] = scipy.stats.zscore(np.log(rips["duration"]))

    ripple_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.POSTPROCESSED_RIPPLES
    )
    rips.to_parquet(ripple_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
