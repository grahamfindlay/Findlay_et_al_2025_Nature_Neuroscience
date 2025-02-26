from ecephys import wne
from ecephys.wne.sglx.pipeline import generate_scoring_bdfs

from findlay2025a import core


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    s3 = core.get_project("shared")
    opts = s3.load_experiment_subject_params(experiment, sglx_subject.name)
    probe = opts["hypnogram_probe"]
    generate_scoring_bdfs.do_experiment_probe(
        experiment,
        probe,
        sglx_subject,
        core.get_project("seahorse"),
        core.get_project("shared"),
    )


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
