from ecephys import wne
from ecephys.wne.sglx.pipeline import generate_probe_sync_table

from findlay2025a import core


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str, stream: str):
    s3 = core.get_project("shared")
    generate_probe_sync_table.do_experiment(experiment, sglx_subject, s3, stream)


def do_project(stream: str):
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment, stream)
