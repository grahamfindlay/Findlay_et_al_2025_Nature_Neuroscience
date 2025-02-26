from ecephys import wne
from ecephys.wne.sglx.pipeline import consolidate_visbrain_hypnograms

from findlay2025a import core


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    for probe in sglx_subject.get_experiment_probes(experiment):
        consolidate_visbrain_hypnograms.do_experiment_probe(
            experiment,
            probe,
            sglx_subject,
            core.get_project("shared"),
            core.get_project("shared"),
        )


def do_project():
    """Process all subjects and experiments in the project."""
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
