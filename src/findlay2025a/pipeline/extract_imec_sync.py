from ecephys import wne
from ecephys.wne.sglx.pipeline import extract_imec_sync

from findlay2025a import core


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str, **kwargs):
    s3 = core.get_project("shared")
    extract_imec_sync.do_experiment(experiment, sglx_subject, s3, **kwargs)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
