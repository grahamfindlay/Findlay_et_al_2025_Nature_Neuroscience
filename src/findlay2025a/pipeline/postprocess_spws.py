import numpy as np
import scipy.stats

from ecephys import wne
from findlay2025a import core, hypnograms, sharp_waves
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    spws = sharp_waves.read_spws(sglx_subject.name, experiment, kind="raw")
    hg = hypnograms.load_consolidated_hypnogram(
        experiment, sglx_subject.name, simplify=True, clean=True
    )
    spws["state"] = hg.get_states(spws["pk_time"])
    spws = spws[spws["duration"] <= 0.350].reset_index(drop=True)
    spws["zlog_amp"] = scipy.stats.zscore(np.log(spws["pk_amp"]))
    spws["zlog_duration"] = scipy.stats.zscore(np.log(spws["duration"]))

    spw_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.POSTPROCESSED_SPWS
    )
    spws.to_parquet(spw_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
