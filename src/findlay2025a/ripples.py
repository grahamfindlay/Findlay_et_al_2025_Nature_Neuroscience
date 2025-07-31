import warnings
from types import MappingProxyType
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ripple_detection.detectors
import xarray as xr
from scipy import signal

from ecephys.npsig.events import detection
from findlay2025a import agg, core
from findlay2025a.constants import Bands, Files

RIPPLE_DATAFRAME_DTYPES = MappingProxyType(
    {
        "start_time": "float64",
        "end_time": "float64",
        "pk_amp": "float64",
        "pk_time": "float64",
        "pk_chan_id": "float64",
        "duration": "float64",
        "pk_chan_y": "float64",
    }
)


def get_empty_ripple_dataframe():
    return pd.DataFrame(columns=RIPPLE_DATAFRAME_DTYPES).astype(RIPPLE_DATAFRAME_DTYPES)


def plot_filter_response(fs: float, w: np.ndarray, h: np.ndarray, title: str):
    "Utility function to plot response functions"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(0.5 * fs * w / np.pi, 20 * np.log10(np.abs(h)))
    ax.set_ylim(-40, 5)
    ax.set_xlim(0, 0.5 * fs)
    ax.grid(True)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Gain (dB)")
    ax.set_title(title)


def ripple_bandpass_filter(
    sampling_frequency: float,
    ripple_band: tuple[float, float] = Bands.RIPPLE,
    plot: bool = False,
) -> tuple[np.ndarray, float]:
    """
    Notes
    ------
    Based on Eric Denovellis' ripple detection package [1].
    [1] https://github.com/Eden-Kramer-Lab/ripple_detection
    """
    ORDER = 101
    nyquist = 0.5 * sampling_frequency
    TRANSITION_BAND = 25
    desired = [
        0,
        ripple_band[0] - TRANSITION_BAND,
        ripple_band[0],
        ripple_band[1],
        ripple_band[1] + TRANSITION_BAND,
        nyquist,
    ]
    taps = signal.remez(ORDER, desired, [0, 1, 0], fs=sampling_frequency)
    if plot:
        w, h = signal.freqz(taps, [1], worN=1024)
        plot_filter_response(sampling_frequency, w, h, "Ripple filter response.")
    return taps, 1.0


def apply_ripple_filter(
    sig: np.ndarray,
    fs: float,
    ripple_band: tuple[float, float] = Bands.RIPPLE,
    plot_filter: bool = False,
) -> np.ndarray:
    """Apply ripple filter to multichannel lfp data.

    Parameters
    ----------
    sig: (n_samples, n_chans)
        The data to filter.
    fs: float
        The sampling frequency of the data.

    Returns
    -------
    filtered_sig: (n_samples, n_chans)
        The filtered signal.

    Notes
    ------
    Based on Eric Denovellis' ripple detection package [1].
    [1] https://github.com/Eden-Kramer-Lab/ripple_detection
    """
    filter_numerator, filter_denominator = ripple_bandpass_filter(
        fs, ripple_band, plot_filter
    )
    is_nan = np.any(np.isnan(sig), axis=-1)
    filtered_sig = np.full_like(sig, np.nan)
    filtered_sig[~is_nan] = signal.filtfilt(
        filter_numerator, filter_denominator, sig[~is_nan], axis=0
    )

    return filtered_sig


def get_ripple_detection_series(
    lfp: xr.DataArray, n_fine_detection_chans: int = 3
) -> xr.DataArray:
    # Filter data in the ripple band
    lff = lfp.copy()
    lff.values = apply_ripple_filter(
        lff.values, lff.fs, Bands.RIPPLE, plot_filter=False
    )  # This takes ~20s

    # Get the timeseries of ripple power that will be used for detection. This takes ~1.5m
    # Instead of a spatial rolling average, we use consensus groups to achieve a similar but better effect
    channels = lff["channel"].values
    n_channels = channels.size
    n_samples = lff["time"].size
    n_consensus_groups = n_channels - n_fine_detection_chans + 1
    consensus_traces = np.zeros((n_samples, n_consensus_groups))
    center_chans = []
    for i in range(n_consensus_groups):
        group_channels = channels[i : i + n_fine_detection_chans]
        center_chans.append(group_channels[len(group_channels) // 2])
        consensus_traces[:, i] = (
            ripple_detection.detectors.get_Kay_ripple_consensus_trace(
                lff.sel({"channel": group_channels}).values, lff.fs
            )
        )
    consensus_traces = xr.DataArray(
        consensus_traces,
        dims=("time", "channel"),
        coords={
            **lff["time"].coords,
            **lff["channel"].sel({"channel": center_chans}).coords,
        },
    )
    return consensus_traces.isel(channel=consensus_traces.argmax(dim="channel"))


# NB: This is the same exact logic as `findlay202a.spws.get_peak_info()`, and
# could be refactored to use that function. I think it's clearer in thise case
# to just duplicate the function, than to create a module just for this one function.
def get_peak_info(sig: xr.DataArray, evts: pd.DataFrame) -> pd.DataFrame:
    """Get properties of each SPW peak.

    Parameters
    ==========
    sig: (time,) DataArray
        The signal to extract peak amplitudes, times, and channels from.
        Probably the series used for event detection.
    evts: DataFrame
        The events, with each event's start and end times.
    """
    evts = evts.copy()
    if len(evts) == 0:
        evts[["pk_amp", "pk_time", "pk_chan_id"]] = np.nan
        return evts

    def _get_peak_info(evt):
        evt_sig = sig.sel(time=slice(evt.start_time, evt.end_time))
        peak = evt_sig[evt_sig.argmax(dim="time")]
        return peak.item(), peak.time.item(), peak.channel.item()

    info = list(map(_get_peak_info, evts.itertuples()))
    evts[["pk_amp", "pk_time", "pk_chan_id"]] = info
    return evts


def detect_ripples(
    detection_lfp,
    n_fine: int,
    detection_threshold: float,
    boundary_threshold: float,
    minimum_duration: float,
    threshold_type: str,
    control_lfp: xr.DataArray = None,
) -> pd.DataFrame:
    if detection_lfp["time"].size < detection_lfp.fs:
        warnings.warn(
            f"Not enough data to apply ripple filter to LFP of shape {detection_lfp.shape}. Skipping ripple detection for this segment."
        )
        return get_empty_ripple_dataframe()

    if threshold_type == "zscore":
        fn = detection.detect_by_zscore
    elif threshold_type == "value":
        fn = detection.detect_by_value
    else:
        raise ValueError(f"thresholdType {threshold_type} not recognized.")

    ser = get_ripple_detection_series(detection_lfp, n_fine_detection_chans=n_fine)
    ripples = fn(
        ser.values,
        ser["time"].values,
        detection_threshold,
        boundary_threshold,
        minimum_duration,
    )

    if control_lfp is not None:
        assert np.all(detection_lfp["time"].values == control_lfp["time"].values), (
            "Control LFP must have the same timestamps as the detection LFP."
        )
        control_ser = get_ripple_detection_series(
            control_lfp, n_fine_detection_chans=n_fine
        )
        control_ripples = detection.detect_by_value(
            control_ser.values,
            control_ser["time"].values,
            ripples.attrs["detection_threshold"],
            ripples.attrs["boundary_threshold"],
            ripples.attrs["minimum_duration"],
        )
        ripples = detection.touch_and_die(ripples, control_ripples)

    ripples.attrs["detection_channel_ids"] = detection_lfp["channel"].values
    ripples.attrs["n_fine"] = n_fine
    ripples = get_peak_info(ser, ripples)
    ripples["duration"] = ripples["end_time"] - ripples["start_time"]
    ripples["pk_chan_y"] = detection_lfp["y"].sel(channel=ripples["pk_chan_id"].values)
    return ripples


def read_ripples(
    subject: str, experiment: str, add_inter_event_intervals=False, kind="postprocessed"
) -> pd.DataFrame:
    nb = core.get_project("seahorse")
    fname = {"raw": Files.RIPPLES, "postprocessed": Files.POSTPROCESSED_RIPPLES}[kind]
    ripple_file = nb.get_experiment_subject_file(experiment, subject, fname)
    ripples = (
        pd.read_parquet(ripple_file).sort_values("start_time").reset_index(drop=True)
    )
    if add_inter_event_intervals:
        ripples["iei"] = ripples.shift(-1)["start_time"] - ripples["end_time"]
    return ripples


def aggregate_ripples() -> Tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    keep = [
        "state",
        "duration",
        "zlog_duration",
        "zlog_amp",
        "freq",
        "zfreq",
        "zdB",
        "zsink",
        "spw_count",
        "spw_zlog_amp",
        "spw_zlog_duration",
        "dB",
    ]
    return agg.aggregate_events(read_ripples, "pk_time", keep)
