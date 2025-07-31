import brainglobe_atlasapi as bgapi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

import ecephys.plot as eplt
import wisc_ecephys_tools as wet
from ecephys.hypnogram import FloatHypnogram
from ecephys.wne.sglx import SGLXSubject
from findlay2025a import agg, core
from wisc_ecephys_tools import rats
from wisc_ecephys_tools.rats.constants import SleepDeprivationExperiments

_ATLAS = bgapi.BrainGlobeAtlas("whs_sd_rat_39um")  # TODO: Cache this properly


class MissingAnatomyError(AttributeError):
    """Exception raised when anatomical data is missing from a recording."""

    pass


def is_cortical(acronym: str) -> bool:
    if acronym == "Cx":
        return True
    if acronym == "HR":
        return False  # Too generic. Could be e.g. PaS.
    elif acronym in _ATLAS.lookup_df.acronym.values:
        ancestors = _ATLAS.get_structure_ancestors(acronym)
        return ("Cx" in ancestors) and ("HF" not in ancestors)
    else:
        return False


def load_cx_ipower(
    experiment: str, subject: str, probe: str, iband: str
) -> xr.DataArray:
    nb = wet.get_sglx_project("shared_nobak")
    f = nb.get_experiment_subject_file(experiment, subject, f"{probe}.{iband}.zarr")
    ipwr = xr.open_dataarray(f)
    # TODO: Instead of using existing acronym coord, it might be better to use
    # existing y/ref_y coords to assign acronym/ref_acronym coords.
    if "acronym" not in ipwr.coords:
        raise MissingAnatomyError(f"No 'acronym' coordinate in {f.name}")
    ipwr = ipwr.isel({"channel": ipwr["acronym"] == ipwr["ref_acronym"]})
    is_cx = [is_cortical(acronym) for acronym in ipwr["acronym"].values]
    if not any(is_cx):
        raise MissingAnatomyError(f"No cortical channels in {f.name}")
    return ipwr.isel({"channel": is_cx})


def _find_first_zero_run(a: np.ndarray, n: int = 1) -> int:
    """Get the index of the first run of n consecutive zeros in a numpy array."""
    first_zero_run = None
    for i in range(len(a) - n + 1):
        if not any(a[i : i + n]):  # All zeros in window
            first_zero_run = i
            break
    return first_zero_run


def _get_threshold_table(
    hist: np.ndarray, bin_edges: np.ndarray, compact: bool = True
) -> pd.DataFrame:
    # Get runs of zeros
    zero_runs = np.split(np.arange(len(hist)), np.where(np.diff(hist == 0))[0] + 1)
    # Filter to only the runs of zeros
    zero_runs = [run for run in zero_runs if hist[run[0]] == 0]
    # Get lengths of zero runs
    run_lengths = [len(run) for run in zero_runs]
    # Get thresholds for different run lengths
    thresholds = {}
    for n in range(1, max(run_lengths) + 1):
        ix = _find_first_zero_run(hist, n)
        if ix is not None:
            thresholds[n] = bin_edges[ix]

    df = pd.DataFrame(
        {"run_length": list(thresholds.keys()), "threshold": list(thresholds.values())}
    )
    if compact:
        # Group by threshold and keep first occurrence (minimum run length)
        df = df.groupby("threshold").first().reset_index()
    return df


def _plot_threshold_table(
    hist: np.ndarray,
    bin_edges: np.ndarray,
    tdf: pd.DataFrame,
    ylim_frac: float = 0.01,
    threshold_ix: int = 0,
) -> tuple[plt.Figure, plt.Axes]:
    fig, axes = plt.subplots(2, 1, figsize=(4, 5), height_ratios=[10, 1])

    axes[0].hist(bin_edges[:-1], bin_edges, weights=hist)
    axes[0].set_ylim(0, ylim_frac * np.max(hist))
    for ix in tdf.index:
        color = "r" if ix == threshold_ix else "k"
        axes[0].axvline(
            tdf.loc[ix, "threshold"], color=color, linestyle="--", linewidth=0.5
        )

    nz = hist > 0
    axes[1].hist(bin_edges[:-1], bin_edges, weights=nz.astype(int))

    return fig, axes


def replace_outliers_kd(
    x: np.ndarray,
    threshold_ix: int = 0,
    bins: int = 1000,
    plot_distribution: bool = False,
    fill_value: float = np.nan,
    **plot_kwargs,
) -> tuple[float, np.ndarray]:
    """
    Usually more conservative (i.e. yields a higher threshold) than even
    np.nanquantile(x, 0.9999).
    """
    hist, bin_edges = np.histogram(x[~np.isnan(x)], bins=bins)
    tdf = _get_threshold_table(hist, bin_edges)
    if plot_distribution:
        fig, axes = _plot_threshold_table(
            hist, bin_edges, tdf, threshold_ix=threshold_ix, **plot_kwargs
        )
    threshold = tdf.loc[threshold_ix, "threshold"]
    x[x > threshold] = fill_value
    return threshold, x


def plot_hypnogram(
    experiment: str,
    subject: SGLXSubject,
    hg: FloatHypnogram,
    state_colors: dict = eplt.publication_colors,
    show_ticklabels: bool = False,
) -> plt.Axes:
    ax = eplt.plot_hypnogram_overlay(
        hg, xlim="hg", figsize=(16, 1), state_colors=state_colors
    )
    wet.rats.cnd_hgs.plot_lights_overlay(
        *wet.rats.cnd_hgs.get_light_dark_periods(experiment, subject),
        ax=ax,
        ymax=1.04,
    )
    if not show_ticklabels:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    return ax


def do_probe(
    experiment: str, subject: str, probe: str, iband: str, plot_ipower: bool = False
) -> xr.DataArray:
    nbsh = wet.get_sglx_project("seahorse")
    _, hg = agg.get_hypnograms(experiment, subject)

    try:
        ipwr = load_cx_ipower(experiment, subject, probe, iband).compute()
    except MissingAnatomyError:
        print(f"No labeled cortex found in {iband} for {experiment} {subject} {probe}")
        return None

    ipwr = ipwr.where(xr.DataArray(hg.covers_time(ipwr["time"]), dims="time"))
    ipwr = ipwr.groupby("acronym").median(dim="channel")
    for acronym in ipwr["acronym"].values:
        da = ipwr.sel({"acronym": acronym})
        _, da.values = replace_outliers_kd(da.values, plot_distribution=True)
        fig = plt.gcf()
        fig.suptitle(acronym)
        fig.savefig(
            nbsh.get_experiment_subject_file(
                experiment, subject, f"{probe}.{acronym}.{iband}_threshold.png"
            )
        )
        plt.close(fig)
        if plot_ipower:
            fs = np.ceil(1 / ipwr["time"].diff(dim="time").median())
            assert fs == 32, f"Sampling rate is {fs} Hz, not 32 Hz"
            sglx_subject = wet.get_sglx_subject(subject)
            ax = plot_hypnogram(experiment, sglx_subject, hg)
            ipwr.sel({"acronym": acronym}).rolling(time=128, center=True).median().plot(
                x="time", color="k", ax=ax
            )
            fig = plt.gcf()
            fig.savefig(
                nbsh.get_experiment_subject_file(
                    experiment,
                    subject,
                    f"{probe}.{acronym}.4s_rolling_median_{iband}.png",
                )
            )
            plt.close(fig)
    savefile = nbsh.get_experiment_subject_file(
        experiment, subject, f"{probe}.cx_{iband}.zarr"
    )
    ipwr.to_zarr(savefile)

    return ipwr


def do_project(iband: str, plot_ipower: bool = False):
    for subject, experiment, probe in rats.utils.get_subject_experiment_probe_tuples(
        experiment_filter=lambda x: x in SleepDeprivationExperiments
    ):
        proj_pairs = [
            (subj, str(expt))
            for subj, expt in core.yield_subject_name_experiment_pairs()
        ]
        if (subject, experiment) not in proj_pairs:
            continue

        print(f"Doing {subject} {experiment} {probe}")
        do_probe(experiment, subject, probe, iband, plot_ipower)
