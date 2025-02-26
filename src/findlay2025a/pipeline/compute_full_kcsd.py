import numpy as np
from ecephys import wne, xrsig

from findlay2025a import core
from findlay2025a.constants import Files


def lazy_mapped_hippocampal_kcsd(subject: str, experiment: str):
    lf = core.open_lfps(subject, experiment)
    nb = core.get_project("seahorse")
    csd_params_file = nb.get_experiment_subject_file(
        experiment, subject, Files.KCSD_PARAMS
    )
    csd_params = np.load(csd_params_file)
    kcsd_kwargs = dict(
        drop=csd_params["channels_omitted"],
        do_lcurve=False,
        gdx=csd_params["gdx"],
        R_init=csd_params["R"],
        lambd=csd_params["lambd"],
    )
    lf = lf.sel(channel=csd_params["channels_used"])

    return xrsig.lazy_mapped_kernel_current_source_density(lf, **kcsd_kwargs)


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    kcsd = lazy_mapped_hippocampal_kcsd(sglx_subject.name, experiment)
    nb = core.get_project("seahorse")
    zarr_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.HIPPOCAMPAL_KCSD
    )
    kcsd.to_zarr(zarr_file, compute=True, mode="w")


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
