from ecephys import wne
from ecephys.wne.sglx.pipeline import emg_from_lfp

from findlay2025a import core


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str, **kwargs):
    emg_from_lfp.do_experiment_probe(
        experiment,
        core.get_spw_probe(experiment, sglx_subject.name),
        sglx_subject,
        core.get_project("shared"),
        core.get_project("shared"),
        core.get_project("seahorse"),
        **kwargs,
    )


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
