from typing import Final

import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr
from ecephys import wne, xrsig

from findlay2025a import core, sharp_waves
from findlay2025a.constants import Experiments, Files

DISABLED_DRIFT: Final[list[tuple[str, Experiments]]] = [
    ("CNPIX6-Eugene", Experiments.NOD),
    ("CNPIX18-Pier", Experiments.COW),
    ("CNPIX19-Otto", Experiments.COW),
]


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    kcsd = core.open_kcsd(sglx_subject.name, experiment)

    # Get SPW detection params
    nb = core.get_project("seahorse")
    spw_params_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.SPW_PARAMS
    )
    spw_params = np.load(spw_params_file)
    spw_detection_kwargs = {
        kw: spw_params[kw]
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
    spw_detection_channels = core.get_shifted_channel_masks(
        shifts, kcsd, spw_params["detection_channel_ids"]
    )
    spw_detection_channels_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.SPW_DETECTION_CHANS
    )
    spw_detection_channels.to_netcdf(spw_detection_channels_file)

    # Get SPWs from estimation period, for their dtypes
    estm_spws_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.ESTM_SPWS
    )
    estm_spws = pd.read_parquet(estm_spws_file)

    # Define computation
    result_dtypes = estm_spws.dtypes.to_dict()

    @dask.delayed
    def _do_block(kcsd: xr.DataArray) -> pd.DataFrame:
        return sharp_waves.detect_sinks(
            kcsd, threshold_type="value", **spw_detection_kwargs
        ).astype(result_dtypes)

    meta = dd.utils.make_meta(result_dtypes)
    results = [
        dd.from_delayed(_do_block(blk_kcsd.sel({"channel": blk_channels})), meta=meta)
        for blk_kcsd, blk_channels in zip(
            xrsig.iterate_timeseries_chunks(kcsd), spw_detection_channels
        )
    ]
    ddf = dd.concat(results)

    # Compute and write to disk
    spw_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.CLASSIC_SPWS
    )
    ddf.to_parquet(spw_file, compute=True)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
