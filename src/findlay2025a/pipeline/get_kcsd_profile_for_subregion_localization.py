import xarray as xr

from ecephys import utils, wne
from findlay2025a import core, hypnograms
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    kcsd_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.HIPPOCAMPAL_KCSD
    )
    kcsd = xr.open_dataarray(kcsd_file, engine="zarr", chunks="auto")

    # Get estimation period
    t1, t2 = core.get_estimation_bounds(sglx_subject, experiment)
    kcsd = kcsd.sel(time=slice(t1, t2))

    # Drop artifacts
    hg = hypnograms.load_consolidated_hypnogram(
        experiment, sglx_subject.name, simplify=True, clean=True
    )
    artifacts = core.load_artifacts(experiment, sglx_subject.name)
    kcsd = core.drop_artifacts(kcsd, artifacts, hg)

    # Compute
    kcsd_profile = xr.Dataset(
        {"minima": kcsd.min(dim="time"), "maxima": kcsd.max(dim="time")}
    )

    # Save
    kcsd_profile_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.ESTM_PERIOD_KCSD_PROFILE
    )
    utils.save_xarray_to_netcdf(kcsd_profile, kcsd_profile_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Processing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
