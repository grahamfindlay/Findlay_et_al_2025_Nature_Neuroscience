import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr


def get_paired_palette():
    PAIRED_PALETTE = sns.color_palette("Paired")
    for i in np.arange(len(PAIRED_PALETTE))[1::2]:
        PAIRED_PALETTE[i] = tuple(
            (np.array(PAIRED_PALETTE[i]) + np.array(PAIRED_PALETTE[i - 1])) / 2
        )
    return PAIRED_PALETTE


def get_nrem_homeostasis_palette():
    p = get_paired_palette()
    return {
        "early_bsl_nrem": p[0],
        "early_rec_nrem_match": p[1],
        "early_ext_wake": p[2],
        "late_ext_wake": p[3],
        "early_rec_nrem": p[4],
        "late_rec_nrem": p[5],
    }


#####
# Plotting traces
#####


def plot_timetrace(da: xr.DataArray, ax: plt.Axes, smoothing: int = 0) -> plt.Axes:
    if smoothing:
        da = da.rolling(time=smoothing).mean()
    sns.lineplot(x=da["time"].values, y=da.values, color="black", linewidth=0.5, ax=ax)
    ax.set(xlabel=None, xmargin=0, xticks=[])


def plot_swa_timetrace(delta: xr.DataArray, ax: plt.Axes, smoothing: int = 0):
    plot_timetrace(delta, ax, smoothing)
    ax.set(ylabel="Cx SWA")


def plot_event_rate_timetrace(events: pd.DataFrame, ax: plt.Axes, **rate_kwargs):
    df = get_smoothed_event_rate(events, **rate_kwargs)
    sns.lineplot(
        x=df["time"].values,
        y=df["rate"].values,
        color="black",
        linewidth=0.5,
        ax=ax,
    )
    ax.set(xlabel=None, ylabel="Rate (Hz)", xmargin=0, xticks=[])


def get_smoothed_event_rate(
    events: pd.DataFrame, time_col: str = "pk_time", smoothing: str = "20s"
) -> pd.DataFrame:
    """Take a dataframe of events, whose times are given by `evt_time_col`, and return the instantanuous event rate at each event time.

    Returns:
    --------
    rate: pd.DataFrame
        Same length as `evts`, with columns `time`, and `rate`, giving the instantanous rate in Hz.
    """
    df = pd.DataFrame(index=pd.to_timedelta(events[time_col], "s")).sort_index()
    df["count"] = 1
    df = df.rolling(smoothing).count()
    df["rate"] = df["count"] / pd.to_timedelta(smoothing).total_seconds()
    df = df.drop(columns="count")
    df = df.reset_index()
    df[time_col] = df[time_col].dt.total_seconds()
    return df.rename(columns={time_col: "time"})


#####
# Plotting histograms over time
#####


def plot_event_amplitude_timehist(
    events: pd.DataFrame,
    ax: plt.Axes,
    binwidth=(180, None),
    time_col="pk_time",
    amp_col="pk_amp",
):
    sns.histplot(
        data=events, x=time_col, y=amp_col, binwidth=binwidth, color="grey", ax=ax
    )
    ax.set(xlabel=None, ylabel="Amplitude", xmargin=0, xticks=[])


def plot_event_duration_timehist(
    events: pd.DataFrame,
    ax: plt.Axes,
    binwidth=(180, 0.010),
    ylim=None,
    time_col="pk_time",
    duration_col="duration",
):
    sns.histplot(
        data=events, x=time_col, y=duration_col, binwidth=binwidth, color="grey", ax=ax
    )
    ax.set(xlabel=None, ylabel="Duration (s)", xmargin=0, xticks=[], ylim=ylim)
