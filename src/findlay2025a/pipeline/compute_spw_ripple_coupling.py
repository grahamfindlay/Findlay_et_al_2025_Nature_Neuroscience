from ecephys import wne

from findlay2025a import core, ripples, sharp_waves
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")

    spws = sharp_waves.read_spws(sglx_subject.name, experiment, kind="postprocessed")
    rips = ripples.read_ripples(sglx_subject.name, experiment, kind="postprocessed")

    spws = core.get_event_coupling(spws, rips, "ripple")
    spw_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.POSTPROCESSED_SPWS
    )
    spws.to_parquet(spw_file)

    rips = core.get_event_coupling(rips, spws, "spw")
    rip_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.POSTPROCESSED_RIPPLES
    )
    rips.to_parquet(rip_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
