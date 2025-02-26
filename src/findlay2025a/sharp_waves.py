from typing import Tuple

import numpy as np
import pandas as pd
import xarray as xr
from ecephys.npsig import event_detection

from findlay2025a import agg, core
from findlay2025a.constants import Files

# -------------------- Detection related functions --------------------


def get_sink_detection_series(
    csd: xr.DataArray, n_fine_detection_chans: int = 3
) -> xr.DataArray:
    """Get the single timeseries to be threshold for current sink detection.

    Parameters
    ==========
    csd: (time, channel) DataArray
        Current source density.
    coarse_detection_chans: DataArray.sel indexer
        Channels used for sink detection.
        Default: Use all channels present in the CSD.
    n_file_detection_chans: int
        The (preferably odd) number of neighboring channels to average over
        when producing CSD estimates for each channel.

    Returns:
    ========
    (time,) DataArray
        The minima, at each time, of the locally smoothed CSD.
        The channel of each minimum is preserved.
    """
    dat = (
        csd.rolling(channel=n_fine_detection_chans, center=True)
        .mean()
        .dropna(dim="channel")
    )
    return -dat.isel(channel=dat.argmin(dim="channel"))


# NB: This is the same exact logic as `findlay2025a.ripples.get_peak_info()`, and
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
        evts[["pk_amp", "pk_time", "pk_chan_id"]] = np.NaN
        return evts

    def _get_peak_info(evt):
        evt_sig = sig.sel(time=slice(evt.start_time, evt.end_time))
        peak = evt_sig[evt_sig.argmax(dim="time")]
        return peak.item(), peak.time.item(), peak.channel.item()

    info = list(map(_get_peak_info, evts.itertuples()))
    evts[["pk_amp", "pk_time", "pk_chan_id"]] = info
    return evts


def detect_sinks(
    csd,
    n_fine,
    detection_threshold,
    boundary_threshold,
    minimum_duration,
    threshold_type,
):
    ser = get_sink_detection_series(csd, n_fine_detection_chans=n_fine)

    if threshold_type == "zscore":
        fn = event_detection.detect_by_zscore
    elif threshold_type == "value":
        fn = event_detection.detect_by_value
    else:
        raise ValueError(f"thresholdType {threshold_type} not recognized.")

    sinks = fn(
        ser.values,
        ser["time"].values,
        detection_threshold,
        boundary_threshold,
        minimum_duration,
    )
    sinks.attrs["detection_channel_ids"] = csd["channel"].values
    sinks.attrs["n_fine"] = n_fine
    sinks = get_peak_info(ser, sinks)
    sinks["duration"] = sinks["end_time"] - sinks["start_time"]
    sinks["pk_chan_y"] = csd["y"].sel(channel=sinks["pk_chan_id"].values)
    return sinks


def read_spws(
    subject: str,
    experiment: str,
    add_inter_event_intervals: bool = False,
    kind: str = "postprocessed",
) -> pd.DataFrame:
    nb = core.get_project("seahorse")
    fname = {"raw": Files.CLASSIC_SPWS, "postprocessed": Files.POSTPROCESSED_SPWS}[kind]
    spw_file = nb.get_experiment_subject_file(experiment, subject, fname)
    spws = pd.read_parquet(spw_file).sort_values("start_time").reset_index(drop=True)
    spws["log_amp"] = np.log10(spws["pk_amp"])
    if add_inter_event_intervals:
        spws["iei"] = spws.shift(-1)["start_time"] - spws["end_time"]
    return spws


def aggregate_spws() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    keep = [
        "state",
        "duration",
        "zlog_duration",
        "pk_amp",
        "log_amp",
        "zlog_amp",
        "ripple_count",
        "ripple_zlog_amp",
        "ripple_zlog_duration",
        "ripple_freq",
        "ripple_zfreq",
        "ripple_zdB",
        "ripple_dB",
    ]
    return agg.aggregate_events(read_spws, "pk_time", keep)
