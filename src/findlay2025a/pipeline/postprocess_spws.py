import numpy as np
import scipy.stats
from ecephys import wne

from findlay2025a import core, sharp_waves
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    s3 = core.get_project("shared")

    hg = s3.load_float_hypnogram(experiment, sglx_subject.name, simplify=True)
    spws = sharp_waves.read_spws(sglx_subject.name, experiment, kind="raw")
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
