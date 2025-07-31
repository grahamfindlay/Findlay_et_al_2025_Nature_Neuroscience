from ecephys import utils, wne
from findlay2025a import core, hypnograms
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    lfp = core.open_lfps(sglx_subject.name, experiment)
    lfp = core.select_hippocampal_channels(experiment, sglx_subject.name, lfp)

    # Get estimation period
    t1, t2 = core.get_estimation_bounds(sglx_subject, experiment)
    lfp = lfp.sel(time=slice(t1, t2))

    # Drop artifacts
    hg = hypnograms.load_consolidated_hypnogram(
        experiment, sglx_subject.name, simplify=True, clean=True
    )
    artifacts = core.load_artifacts(experiment, sglx_subject.name)
    lfp = core.drop_artifacts(lfp, artifacts, hg)

    # Compute
    relative_ripple_power = core.get_relative_ripple_power_profile(lfp)

    # Save
    ripple_profile_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.ESTM_PERIOD_RIPPLE_PROFILE
    )
    utils.save_xarray_to_netcdf(relative_ripple_power, ripple_profile_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
