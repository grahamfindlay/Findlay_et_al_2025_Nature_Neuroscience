import shutil

import dask
import dask.dataframe as dd
import matplotlib.pyplot as plt
import msalign
import numpy as np
import pandas as pd
import scipy.signal
import xarray as xr
from ecephys import wne, xrsig

from findlay2025a import core
from findlay2025a.constants import Files

PROFILE_SMOOTHING = (
    5  # Alignment is greatly improved when the power profile is smoothed first
)
SHIFT_RANGE = [-10, 10]  # Maximum and minimum channel shifts to check
MSALIGN_WIDTH = 10  # Width of the gaussian peak used for alignment, in channels


def get_channel_shift_using_msalign(
    tst: xr.DataArray, ref: xr.DataArray, prominence_pct=0.05, y=None, plot=False
) -> int:
    """Align the test profile to the reference profile by shifting channels.
    Uses peak finding and alignment, with alignment priority weighted by peak prominence.
    The shift value return is the amount that the reference profile should be shifted
    using np.roll / scipy.shift to match the test profile (So the opposite of the
    shift value is the amount that the test profile must be shifted to match the reference).
    Promimence PCT can be used so that only larger peaks are used for alignment.
    """
    # Get peak info for alignment
    if prominence_pct is not None:

        def thresh(x):
            return float(np.abs(np.max(x) - np.min(x)) * prominence_pct)

        pks = scipy.signal.find_peaks(ref, prominence=thresh(ref))[0]
    else:
        pks = scipy.signal.find_peaks(ref)[0]
    proms = scipy.signal.peak_prominences(ref, pks)[0]
    ixs = np.arange(ref.size)
    # Actually align
    aligned, shift = msalign.msalign(
        ixs,
        np.atleast_2d(tst),
        pks,
        weights=proms,
        width=MSALIGN_WIDTH,  # Width of the gaussian peak used for alignment, in channels
        shift_range=SHIFT_RANGE,
        only_shift=True,
        return_shifts=True,
        align_by_index=True,
    )

    if plot:
        # Plot reference and test signals, and the realigned test signal.
        aligned[aligned == 0] = np.nan  # Just for plotting, so y axis doesn't go to 0
        fig, ax = plt.subplots(figsize=(32, 6))
        if y is None:
            y = ixs
        ax.plot(y, ref)
        ax.plot(y, tst, c="r", alpha=0.5)
        ax.plot(y, aligned.squeeze(), c="g")
        for pk in pks:
            ax.axvline(y[pk], c="k", linestyle=":")

    return shift.squeeze()


def do_experiment(
    sglx_subject: wne.sglx.SGLXSubject, experiment: str, estimation_chunk_size=None
):
    """
    estimation_chunk_size: int
        In seconds. We will get 1 shift estimate per chunk, so this is your resolution.
        For ~140s chunks, takes ~4.5m on 48h of data, and works well.
        For 900s chunks, takes ~2.5m on 48h of data, and works well.
    """
    # Get the ripple power profile from the estimation period
    # All shift values will be relative to this period
    nb = core.get_project("seahorse")
    ripple_profile_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.ESTM_PERIOD_RIPPLE_PROFILE
    )
    estm_period_profile = xr.open_dataarray(ripple_profile_file)

    # Open LFPs, and rechunk. We get 1 shift estimate per chunk, so chunksize determines resoluton
    lfp = core.open_hippocampal_lfps(sglx_subject.name, experiment)
    if estimation_chunk_size is not None:
        lfp = lfp.chunk(chunks={"time": int(lfp.fs * estimation_chunk_size)})

    # Get shift for each chunk/block.
    result_dtypes = {
        "block_shift": int,
        "block_start_time": float,
        "block_end_time": float,
    }

    @dask.delayed
    def _get_block_shift_info(blk_lfp: xr.DataArray, reference_profile: xr.DataArray):
        block_profile = core.get_relative_ripple_power_profile(blk_lfp)
        block_shift = get_channel_shift_using_msalign(
            block_profile.rolling(channel=PROFILE_SMOOTHING).mean().dropna("channel"),
            reference_profile.rolling(channel=PROFILE_SMOOTHING)
            .mean()
            .dropna("channel"),
        )
        return pd.DataFrame(
            {
                "block_shift": int(block_shift),
                "block_start_time": float(blk_lfp["time"].values[0]),
                "block_end_time": float(blk_lfp["time"].values[-1]),
            },
            index=[0],
        ).astype(result_dtypes)

    meta = dd.utils.make_meta(result_dtypes)
    results = [
        dd.from_delayed(_get_block_shift_info(blk_lfp, estm_period_profile), meta=meta)
        for blk_lfp in xrsig.iterate_timeseries_chunks(lfp)
    ]
    ddf = dd.concat(results)

    # Compute and write to disk
    shift_info_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_POWER_SHIFT
    )
    if shift_info_file.exists():
        shutil.rmtree(shift_info_file)
    ddf.to_parquet(shift_info_file, compute=True)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Processing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
