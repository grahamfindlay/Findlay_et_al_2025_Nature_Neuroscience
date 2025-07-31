from types import MappingProxyType
from typing import Final

import numpy as np

from ecephys import wne
from findlay2025a import core, hypnograms, ripples
from findlay2025a.constants import Files

RIPPLE_DETECTION_RADIUS: Final = (
    80  # Area round the ripple peak channel to search, in microns
)

RIPPLE_DETECTION_PARAMS: Final = MappingProxyType(
    dict(
        n_fine=3,
        detection_threshold=2.0,
        boundary_threshold=0.5,
        minimum_duration=0.020,
        threshold_type="zscore",
    )
)


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment):
    # Get LFP
    t1, t2 = core.get_estimation_bounds(sglx_subject, experiment)
    lfp = core.open_hippocampal_lfps(sglx_subject.name, experiment).sel(
        time=slice(t1, t2)
    )
    wm_lfp = core.open_white_matter_lfps(
        sglx_subject.name, experiment, drop_bad_channels=True
    ).sel(time=slice(t1, t2))

    # Drop artifacts
    hg = hypnograms.load_consolidated_hypnogram(
        experiment, sglx_subject.name, simplify=True, clean=True
    )
    artifacts = core.load_artifacts(experiment, sglx_subject.name)
    lfp = core.drop_artifacts(lfp, artifacts, hg)
    wm_lfp = core.drop_artifacts(wm_lfp, artifacts, hg)

    # Get Ripple thresholds
    _, ripple_detection_channels = core.get_hippocampal_subregion_channels(
        sglx_subject.name,
        experiment,
        lfp,
        "ripple_power_peak",
        radius=RIPPLE_DETECTION_RADIUS,
    )
    detection_lfp = lfp.sel(channel=ripple_detection_channels).compute()
    control_lfp = wm_lfp.compute()
    rips = ripples.detect_ripples(
        detection_lfp, **RIPPLE_DETECTION_PARAMS, control_lfp=control_lfp
    )
    nb = core.get_project("seahorse")
    ripple_params_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_PARAMS
    )
    np.savez(ripple_params_file, **rips.attrs)
    estm_ripples_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.ESTM_RIPPLES
    )
    rips.to_parquet(estm_ripples_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
