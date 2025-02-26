from ecephys import wne
from ecephys.wne.sglx.pipeline import get_scoring_signals

from findlay2025a import core


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    get_scoring_signals.do_experiment(
        experiment,
        sglx_subject,
        core.get_project("shared"),
        core.get_project("seahorse"),
    )


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
