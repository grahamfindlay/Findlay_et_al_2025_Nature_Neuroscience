from types import MappingProxyType
from typing import Final

import numpy as np
from ecephys import wne

from findlay2025a import core, sharp_waves
from findlay2025a.constants import Files

SPW_DETECTION_RADIUS: Final[float] = 100  # in Microns

SPW_ESTIMATION_PARAMS: Final[MappingProxyType] = MappingProxyType(
    dict(
        n_fine=3,
        detection_threshold=2,
        boundary_threshold=1.5,
        minimum_duration=0.020,
        threshold_type="zscore",
    )
)


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    s3 = core.get_project("shared")

    t1, t2 = core.get_estimation_bounds(sglx_subject, experiment)
    kcsd = core.open_kcsd(sglx_subject.name, experiment).sel(time=slice(t1, t2))

    # Drop artifacts
    hg = s3.load_float_hypnogram(experiment, sglx_subject.name, simplify=True)
    artifacts = core.load_artifacts(experiment, sglx_subject.name)
    kcsd = core.drop_artifacts(kcsd, artifacts, hg)

    _, spw_detection_channels = core.get_hippocampal_subregion_channels(
        sglx_subject.name,
        experiment,
        kcsd,
        "stratum_radiatum_peak",
        radius=SPW_DETECTION_RADIUS,
    )
    spws = sharp_waves.detect_sinks(
        kcsd.sel(channel=spw_detection_channels).compute(), **SPW_ESTIMATION_PARAMS
    )
    spw_params_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.SPW_PARAMS
    )
    np.savez(spw_params_file, **spws.attrs)
    estm_spws_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.ESTM_SPWS
    )
    spws.to_parquet(estm_spws_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
