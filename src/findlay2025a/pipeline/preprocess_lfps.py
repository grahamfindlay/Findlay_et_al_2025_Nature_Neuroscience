from ecephys import wne
from ecephys.wne.sglx.pipeline import preprocess_lfps

from findlay2025a import core


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    preprocess_lfps.do_experiment(
        experiment,
        sglx_subject,
        core.get_project("shared"),
        core.get_project("seahorse"),
        core.get_project("shared"),
    )


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
