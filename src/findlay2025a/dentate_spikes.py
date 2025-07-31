from typing import Final, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal
import scipy.stats
import xarray as xr

from ecephys import dasig, xrsig
from findlay2025a import agg, core
from findlay2025a.constants import Files

DETECTION_LOWCUT: Final = 5
DETECTION_HIGHCUT: Final = 100
DETECTION_FILTER_ORDER: Final = 3


def extract_dentate_spikes(
    df: pd.DataFrame,
    reference_spikes: pd.DataFrame = None,
    dvorak_width: tuple[float, float] = (0.005, 0.015),
    jones_width: tuple[float, float] = (0.000, np.inf),
    dvorak_zlog_prominence: float = 1.5,
    findlay_zlog_prominence: float = 1.5,
) -> pd.DataFrame:
    if reference_spikes is None:
        satisfies_prominence = (
            df["dvorak_zlog_prominence"] > dvorak_zlog_prominence
        ) & (df["findlay_zlog_prominence"] > findlay_zlog_prominence)
    else:
        satisfies_prominence = (
            df["dvorak_prominence"] >= reference_spikes["dvorak_prominence"].min()
        ) & (df["findlay_prominence"] >= reference_spikes["findlay_prominence"].min())
    satisfies_width = (df["dvorak_width"] >= dvorak_width[0]) & (
        df["dvorak_width"] <= dvorak_width[1]
    )
    satisfies_width &= (df["jones_width"] >= jones_width[0]) & (
        df["jones_width"] <= jones_width[1]
    )

    return df.loc[satisfies_prominence & satisfies_width].reset_index(drop=True)


def read_dspks(
    subject: str, experiment: str, kind="postprocessed", **raw_kwargs
) -> pd.DataFrame:
    nb = core.get_project("seahorse")
    if kind == "raw":
        dpk_file = nb.get_experiment_subject_file(experiment, subject, Files.DPKS)
        ddf = pd.read_parquet(dpk_file).sort_values("peak_time").reset_index(drop=True)
        dspks = extract_dentate_spikes(ddf, **raw_kwargs)
    elif kind == "postprocessed":
        dspk_file = nb.get_experiment_subject_file(
            experiment, subject, Files.POSTPROCESSED_DSPKS
        )
        dspks = (
            pd.read_parquet(dspk_file).sort_values("peak_time").reset_index(drop=True)
        )
    else:
        raise ValueError(f"Unknown kind: {kind}")
    return dspks


def aggregate_dspks() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    keep = [
        "state",
        "zlog_height",
        "dvorak_zlog_prominence",
        "findlay_zlog_prominence",
    ]
    return agg.aggregate_events(read_dspks, "peak_time", keep)


def get_region_bounds(da: xr.DataArray) -> pd.DataFrame:
    channel_info = da.channel.to_dataframe()
    return channel_info.groupby("acronym")["y"].describe()[["min", "max"]]


def get_leaf_detection_series(
    sink_zone_csd: xr.DataArray, source_zone_csd: xr.DataArray
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    sink_series = sink_zone_csd.argmin(dim="channel")
    sink_series = sink_zone_csd.isel(channel=sink_series)

    source_series = source_zone_csd.argmax(dim="channel")
    source_series = source_zone_csd.isel(channel=source_series)

    detection_series = -(sink_series - source_series)

    return sink_series, source_series, detection_series


def detect_peaks(
    detection_series: xr.DataArray,
    fs: float,
    dorsal_sink_series: xr.DataArray = None,
    dorsal_source_series: xr.DataArray = None,
    ventral_sink_series: xr.DataArray = None,
    ventral_source_series: xr.DataArray = None,
) -> pd.DataFrame:
    pks, pk_props = scipy.signal.find_peaks(detection_series)  # Get peaks
    trs, tr_props = scipy.signal.find_peaks(-detection_series)  # Get troughs

    # Keep only peaks with preceeding and succeeding troughs
    keep = (pks > trs[0]) & (pks < trs[-1])
    pks = pks[keep]

    # Get the troughs that preceed and succeed each peak
    i = np.searchsorted(trs, pks, side="left")
    assert np.all(i > 0), "Found a peak without a preceeding trough"
    assert np.all(i < trs.size), "Found a peak without a succeeding trough"
    left_bases = trs[i - 1]  # For each peak, get its left base
    right_bases = trs[i]  # For each peak, get its right base

    # For each peak, determine which of its bases is nearest in time
    rise_durations = pks - left_bases
    fall_durations = right_bases - pks
    nearest_bases = np.where(
        np.argmin(np.stack([rise_durations, fall_durations]), axis=0),
        right_bases,
        left_bases,
    )

    # For each peak, determine which of its bases is highest and lowest
    left_base_amps = detection_series.values[left_bases]
    right_base_amps = detection_series.values[right_bases]
    deepest_bases = np.where(
        np.argmin(np.stack([left_base_amps, right_base_amps]), axis=0),
        right_bases,
        left_bases,
    )
    highest_bases = np.where(
        np.argmax(np.stack([left_base_amps, right_base_amps]), axis=0),
        right_bases,
        left_bases,
    )

    # For each peak, get its width and prominence as defined by Jones et al. 2021
    pk_amps = detection_series.values[pks]
    jones_proms = pk_amps - detection_series.values[nearest_bases]
    jones_width = right_bases - left_bases

    # For each peak, get its width and prominence as defined by Dvorak et al. 2021
    dvorak_proms = pk_amps - detection_series.values[deepest_bases]
    dvorak_width = np.abs(pks - nearest_bases) * 2

    # For each peak, get its width and prominence using the highest base.
    findlay_proms = pk_amps - detection_series.values[highest_bases]

    # For each peak, get its width and prominence using the scipy method
    scipy_proms, scipy_left_bases, scipy_right_bases = scipy.signal.peak_prominences(
        detection_series.values, pks
    )
    scipy_width = scipy_right_bases - scipy_left_bases

    df = pd.DataFrame(
        {
            "peak_frame": pks,
            "peak_time": detection_series.time.values[pks],
            "peak_height": pk_amps,
            "dvorak_prominence": dvorak_proms,
            "dvorak_width": dvorak_width / fs,
            "jones_prominence": jones_proms,
            "jones_width": jones_width / fs,
            "findlay_prominence": findlay_proms,
            "scipy_prominence": scipy_proms,
            "scipy_width": scipy_width / fs,
        }
    )
    if dorsal_sink_series is not None:
        df["dorsal_sink_channel"] = dorsal_sink_series.isel(time=pks).channel.values
        df["dorsal_sink_y"] = dorsal_sink_series.isel(time=pks).y.values
        df["dorsal_sink_structure"] = dorsal_sink_series.isel(time=pks).acronym.values
    if dorsal_source_series is not None:
        df["dorsal_source_channel"] = dorsal_source_series.isel(time=pks).channel.values
        df["dorsal_source_y"] = dorsal_source_series.isel(time=pks).y.values
    if ventral_sink_series is not None:
        df["ventral_sink_channel"] = ventral_sink_series.isel(time=pks).channel.values
        df["ventral_sink_y"] = ventral_sink_series.isel(time=pks).y.values
        df["ventral_sink_structure"] = ventral_sink_series.isel(time=pks).acronym.values
    if ventral_source_series is not None:
        df["ventral_source_channel"] = ventral_source_series.isel(
            time=pks
        ).channel.values
        df["ventral_source_y"] = ventral_source_series.isel(time=pks).y.values

    if dorsal_sink_series is not None and dorsal_source_series is not None:
        df["dorsal_sink_source_separation"] = (
            df["dorsal_sink_y"] - df["dorsal_source_y"]
        )
        df = df.loc[df["dorsal_sink_source_separation"] > 0.0]
    if ventral_sink_series is not None and ventral_source_series is not None:
        df["ventral_sink_source_separation"] = (
            df["ventral_source_y"] - df["ventral_sink_y"]
        )
        df = df.loc[df["ventral_sink_source_separation"] > 0.0]

    if dorsal_sink_series is not None and ventral_sink_series is not None:
        df["sink_separation"] = df["dorsal_sink_y"] - df["ventral_sink_y"]
    if dorsal_source_series is not None and ventral_source_series is not None:
        df["source_separation"] = df["dorsal_source_y"] - df["ventral_source_y"]

    df = add_zlog_prominence(df)
    return df.reset_index(drop=True)


def add_zlog_prominence(df: pd.DataFrame) -> pd.DataFrame:
    df["dvorak_zlog_prominence"] = scipy.stats.zscore(np.log10(df["dvorak_prominence"]))
    df["jones_zlog_prominence"] = scipy.stats.zscore(np.log10(df["jones_prominence"]))
    df["scipy_zlog_prominence"] = scipy.stats.zscore(np.log10(df["scipy_prominence"]))
    df["findlay_zlog_prominence"] = scipy.stats.zscore(
        np.log10(df["findlay_prominence"])
    )
    return df


def examine_event(
    evts: pd.DataFrame,
    lfp: xr.DataArray,
    kcsd: xr.DataArray,
    detection_series: xr.DataArray = None,
    i: int = None,
):
    if i is None:
        i = np.random.choice(len(evts))
    evt = evts.iloc[i]

    start = evt.peak_time - 1.0
    end = evt.peak_time + 1.0

    if detection_series is not None:
        fig, axes = plt.subplots(
            3, 1, figsize=(32, 15), sharex=True, height_ratios=[6, 6, 1]
        )
    else:
        fig, axes = plt.subplots(2, 1, figsize=(32, 14), sharex=True)

    kcsd.sel(time=slice(start, end)).plot(
        x="time", y="y", add_colorbar=False, ax=axes[0]
    )
    xrsig.plot_traces(
        lfp.sel(time=slice(start, end)), ax=axes[1], vspace=200, flip_dv=True
    )
    axes[1].axvline(evt.peak_time, color="r", ls=":")

    if detection_series is not None:
        detection_series.sel(time=slice(start, end)).plot(x="time", ax=axes[2])
        detection_series.isel(time=evts.peak_frame.values).sel(
            time=slice(start, end)
        ).plot.scatter(x="time", ax=axes[2], color="r")

    return evt


def shift_blocks(da: xr.DataArray, shifts: pd.DataFrame) -> xr.DataArray:
    xrsig.validate_2d_timeseries(da)
    block_shifts = core.get_block_shifts(da, shifts)
    res = da.copy()
    res.data = dasig.shift_blocks(
        res.data, block_shifts, axis=res.get_axis_num("channel")
    )
    return res


def visualize_detection_zones(
    start: float,
    end: float,
    kcsd: xr.DataArray,
    lfp: xr.DataArray,
    region_bounds: pd.DataFrame,
    dorsal_sink_chans: list[int],
    dorsal_source_chans: list[int],
    ventral_sink_chans: list[int],
    ventral_source_chans: list[int],
):
    ls1 = (0, (1, 10))
    ls2 = (0, (1, 5))
    fig, axes = plt.subplots(2, 1, figsize=(32, 14), sharex=True)
    kcsd.sel(time=slice(start, end)).plot(
        x="time", y="y", add_colorbar=False, ax=axes[0]
    )
    for region in region_bounds.itertuples():
        axes[0].axhline(region.min, c="k", linestyle=ls1)
        axes[0].axhline(region.max, c="k", linestyle=ls1)

    if dorsal_sink_chans is not None:
        dorsal_sink_lo = float(kcsd.sel(channel=dorsal_sink_chans)["y"].min())
        dorsal_sink_hi = float(kcsd.sel(channel=dorsal_sink_chans)["y"].max())
        axes[0].axhline(dorsal_sink_lo, c="g", linestyle=ls2)
        axes[0].axhline(dorsal_sink_hi, c="g", linestyle=ls2)

    if dorsal_source_chans is not None:
        dorsal_source_lo = float(kcsd.sel(channel=dorsal_source_chans)["y"].min())
        dorsal_source_hi = float(kcsd.sel(channel=dorsal_source_chans)["y"].max())
        axes[0].axhline(dorsal_source_lo, c="purple", linestyle=ls2)
        axes[0].axhline(dorsal_source_hi, c="purple", linestyle=ls2)

    if ventral_sink_chans is not None:
        ventral_sink_lo = float(kcsd.sel(channel=ventral_sink_chans)["y"].min())
        ventral_sink_hi = float(kcsd.sel(channel=ventral_sink_chans)["y"].max())
        axes[0].axhline(ventral_sink_lo, c="g", linestyle=ls2)
        axes[0].axhline(ventral_sink_hi, c="g", linestyle=ls2)

    if ventral_source_chans is not None:
        ventral_source_lo = float(kcsd.sel(channel=ventral_source_chans)["y"].min())
        ventral_source_hi = float(kcsd.sel(channel=ventral_source_chans)["y"].max())
        axes[0].axhline(ventral_source_lo, c="purple", linestyle=ls2)
        axes[0].axhline(ventral_source_hi, c="purple", linestyle=ls2)

    xrsig.plot_traces(
        lfp.sel(time=slice(start, end)),
        ax=axes[1],
        vspace=200,
        flip_dv=True,
        chan_labels="acronym",
        chan_colors="acronym",
    )
