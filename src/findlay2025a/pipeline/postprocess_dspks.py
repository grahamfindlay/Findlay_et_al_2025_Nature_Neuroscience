import numpy as np
import scipy.stats

from ecephys import wne
from findlay2025a import core, dentate_spikes, hypnograms
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    hg = hypnograms.load_consolidated_hypnogram(
        experiment, sglx_subject.name, simplify=True, clean=True
    )
    dspks = dentate_spikes.read_dspks(sglx_subject.name, experiment, kind="raw")
    dspks["state"] = hg.get_states(dspks["peak_time"])
    dspks["zlog_height"] = scipy.stats.zscore(np.log(dspks["peak_height"]))
    dspk_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.POSTPROCESSED_DSPKS
    )
    dspks.to_parquet(dspk_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
