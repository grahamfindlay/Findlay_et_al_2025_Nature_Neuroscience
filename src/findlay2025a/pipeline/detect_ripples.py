from typing import Final

import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr
from ecephys import wne, xrsig

from findlay2025a import core, ripples
from findlay2025a.constants import Experiments, Files

DISABLED_DRIFT: Final[list[tuple[str, Experiments]]] = [
    ("CNPIX6-Eugene", Experiments.NOD),
    ("CNPIX18-Pier", Experiments.COW),
    ("CNPIX19-Otto", Experiments.COW),
]


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    lfp = core.open_hippocampal_lfps(sglx_subject.name, experiment)
    wm_lfp = core.open_white_matter_lfps(
        sglx_subject.name, experiment, drop_bad_channels=True
    )

    # Get ripple detection params
    nb = core.get_project("seahorse")
    ripple_params_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_PARAMS
    )
    ripple_params = np.load(ripple_params_file)
    ripple_detection_kwargs = {
        kw: ripple_params[kw]
        for kw in [
            "detection_threshold",
            "boundary_threshold",
            "minimum_duration",
            "n_fine",
        ]
    }
    shifts = core.read_shifts(sglx_subject.name, experiment)
    if (sglx_subject.name, experiment) in DISABLED_DRIFT:
        shifts["block_shift"] = 0
    ripple_detection_channels = core.get_shifted_channel_masks(
        shifts, lfp, ripple_params["detection_channel_ids"]
    )
    ripple_detection_channels_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_DETECTION_CHANS
    )
    ripple_detection_channels.to_netcdf(ripple_detection_channels_file)

    # Get ripples from estimation period, for their dtypes
    estm_ripples_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.ESTM_RIPPLES
    )
    estm_ripples = pd.read_parquet(estm_ripples_file)

    # Define computation
    result_dtypes = estm_ripples.dtypes.to_dict()

    @dask.delayed
    def _do_block(lfp: xr.DataArray, ctrl_lfp: xr.DataArray) -> pd.DataFrame:
        return ripples.detect_ripples(
            lfp, threshold_type="value", **ripple_detection_kwargs, control_lfp=ctrl_lfp
        ).astype(result_dtypes)

    meta = dd.utils.make_meta(result_dtypes)
    results = [
        dd.from_delayed(
            _do_block(blk_lfp.sel({"channel": blk_channels}), blk_ctrl_lfp), meta=meta
        )
        for blk_lfp, blk_ctrl_lfp, blk_channels in zip(
            xrsig.iterate_timeseries_chunks(lfp),
            xrsig.iterate_timeseries_chunks(wm_lfp),
            ripple_detection_channels,
        )
    ]
    ddf = dd.concat(results)

    # Compute and write to disk
    ripples_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLES
    )
    ddf.to_parquet(ripples_file, compute=True)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
