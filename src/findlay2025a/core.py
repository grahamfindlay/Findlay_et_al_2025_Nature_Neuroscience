import pathlib
import warnings
from types import MappingProxyType
from typing import Final, Generator, List, Tuple

import numpy as np
import pandas as pd
import scipy.signal
import xarray as xr

import wisc_ecephys_tools as wet
from ecephys import hypnogram, utils, wne, xrsig
from ecephys.wne.constants import Files as WNEFiles
from findlay2025a.constants import Bands, Files
from findlay2025a.constants import Experiments as Exps

_WNE_PROJECTS = {
    "seahorse": wet.get_sglx_project("seahorse"),
    "shared": wet.get_sglx_project("shared"),
}

MANIFEST = MappingProxyType(
    {
        "CNPIX2-Segundo": [Exps.NOD],
        "CNPIX3-Valentino": [Exps.NOD],
        "CNPIX4-Doppio": [Exps.NOD],
        "CNPIX5-Alessandro": [Exps.NOD],
        "CNPIX6-Eugene": [Exps.NOD],
        "CNPIX8-Allan": [Exps.NOD],
        "CNPIX9-Luigi": [Exps.NOD],
        "CNPIX10-Charles": [Exps.NOD],
        "CNPIX11-Adrian": [Exps.NOD, Exps.COW],
        "CNPIX12-Santiago": [Exps.NOD, Exps.COW],
        "CNPIX14-Francis": [Exps.NOD, Exps.COW],
        "CNPIX15-Claude": [Exps.NOD, Exps.COW, Exps.CTN],
        "CNPIX17-Hans": [Exps.NOD, Exps.COW, Exps.CTN],
        "CNPIX18-Pier": [Exps.NOD, Exps.COW, Exps.CTN],
        "CNPIX19-Otto": [Exps.NOD, Exps.COW, Exps.CTN],
        "CNPIX20-Ernst": [Exps.NOD, Exps.COW, Exps.CTN],
    }
)

SIMPLIFIED_STATES = MappingProxyType(
    {
        "Wake": "Wake",
        "W": "Wake",
        "aWk": "Wake",
        "qWk": "Wake",
        "QWK": "Wake",
        "Arousal": "MA",
        "MA": "MA",
        "Trans": "Other",
        "NREM": "NREM",
        "N1": "NREM",
        "N2": "NREM",
        "IS": "IS",
        "REM": "REM",
        "Art": "Artifact",
        "Artifact": "Artifact",
        "scrambling": "Artifact",
        "None": "Other",
    }
)
ARTIFACT_STATES: Final[List[str]] = [
    k for k, v in SIMPLIFIED_STATES.items() if v == "Artifact"
]

EXPERIMENT_DISPLAY_NAMES = {
    Exps.NOD: "Novelty",
    Exps.COW: "Locomotion",
    Exps.CTN: "Dual",
}

CONDITION_DISPLAY_NAMES = {
    "early_bsl_nrem": "Early.BSL.NREM",
    "early_rec_nrem_match": "Late.BSL.NREM",
    "early_ext_wake": "Early.EXT.Wake",
    "late_ext_wake": "Late.EXT.Wake",
    "early_rec_nrem": "Early.REC.NREM",
    "late_rec_nrem": "Late.REC.NREM",
    "early_nod_wake": "Early.NOD.Wake",
    "early_sd_wake": "Early.SD.Wake",
    "late_sd_wake": "Late.SD.Wake",
    "sd_wake": "SD.Wake",
    "ext_wake": "EXT.Wake",
    "nod_wake": "NOD.Wake",
    "cow_wake": "COW.Wake",
    "ctn_wake": "CTN.Wake",
    "late_cow_wake": "Late.COW.Wake",
}


# TODO: Replace with wet.get_sglx_project()
def get_project(name: str = "seahorse") -> wne.sglx.SGLXProject:
    return _WNE_PROJECTS[name]


def get_subjects(experiment: str) -> list[str]:
    return [sub for sub, exps in MANIFEST.items() if experiment in exps]


def yield_subject_name_experiment_pairs(
    experiments=[Exps.NOD, Exps.COW, Exps.CTN],
) -> Generator[tuple[str, str], None, None]:
    for experiment in experiments:
        for subject in get_subjects(experiment):
            yield (subject, experiment)


def yield_sglx_subject_experiment_pairs(
    experiments=[Exps.NOD, Exps.COW, Exps.CTN],
) -> Generator[tuple[wne.sglx.SGLXSubject, str], None, None]:
    for experiment in experiments:
        for subject in get_subjects(experiment):
            yield (wet.get_sglx_subject(subject), experiment)


# TODO: Replace with Project.load_experiment_subject_params()
def get_experiment_params(experiment: str, subject: str) -> dict:
    s3 = get_project("shared")
    return s3.load_experiment_subject_params(experiment, subject)


def get_spw_probe(experiment: str, subject: str) -> str:
    opts = get_experiment_params(experiment, subject)
    return opts["spwrProbe"]


def get_estimation_bounds(
    wne_subject: wne.sglx.SGLXSubject, experiment: str
) -> tuple[float, float]:
    params = get_experiment_params(experiment, wne_subject.name)
    start_datetime = pd.to_datetime(params["estmPeriod"][0])
    end_datetime = pd.to_datetime(params["estmPeriod"][1])
    start_time = wne_subject.dt2t(experiment, "imec0", start_datetime)
    end_time = wne_subject.dt2t(experiment, "imec0", end_datetime)
    # Using imec0 here is a bit of a hack, to ensure that the times are in the
    # canonical timebase. It should not matter which probe is used.
    # The biggest drawback of this approach is that imec0 needs to be available.
    # If it were not, we'd have to:
    # 1. Take `probe` as an argument.
    # 2. Create a time synchronizer object here and use that.
    return start_time, end_time


#####
# Functions for opening, loading, and selecting data
#####


def open_kcsd(subject: str, experiment: str, chunks="auto") -> xr.DataArray:
    nb = get_project("seahorse")
    kcsd_file = nb.get_experiment_subject_file(
        experiment, subject, Files.HIPPOCAMPAL_KCSD
    )
    kcsd = xr.open_dataarray(kcsd_file, engine="zarr", chunks=chunks)
    kcsd.attrs["fs"] = None  # Just used to pass xrsig validation
    return kcsd


def open_lfps(
    subject: str, experiment: str, drop_bad_channels: bool = False, **kwargs
) -> xr.DataArray:
    s3 = get_project("shared")
    probe = get_spw_probe(experiment, subject)
    badchan_proj = s3 if drop_bad_channels else None
    return wne.utils.open_lfps(
        get_project("seahorse"),
        subject,
        experiment,
        probe,
        anatomy_proj=s3,
        badchan_proj=badchan_proj,
        **kwargs,
    )


def get_hippocampal_bounds(experiment: str, subject: str):
    s3 = get_project("shared")
    params = s3.load_experiment_subject_params(experiment, subject)
    probe = params["spwrProbe"]
    [lo, hi] = params["probes"][probe]["structureBounds"]["hippocampus"]
    return (lo, hi)


def select_hippocampal_channels(
    experiment: str, subject: str, da: xr.DataArray
) -> xr.DataArray:
    lo, hi = get_hippocampal_bounds(experiment, subject)
    is_hippocampal = (da["y"] >= lo) & (da["y"] <= hi)
    return da.sel(channel=is_hippocampal)


def open_hippocampal_lfps(subject: str, experiment: str, **kwargs) -> xr.DataArray:
    lfps = open_lfps(subject, experiment, **kwargs)
    return select_hippocampal_channels(experiment, subject, lfps)


def select_cortical_channels(
    experiment: str,
    subject: str,
    da: xr.DataArray,
    wm_thickness: float = 200,
) -> xr.DataArray:
    s3 = get_project("shared")
    params = s3.load_experiment_subject_params(experiment, subject)
    probe = params["spwrProbe"]
    [hc_lo, hc_hi] = params["probes"][probe]["structureBounds"]["hippocampus"]
    is_cortical = da["y"] > (hc_hi + wm_thickness)
    return da.sel({"channel": is_cortical})


def open_cortical_lfps(subject: str, experiment: str, **kwargs) -> xr.DataArray:
    lfps = open_lfps(subject, experiment, **kwargs)
    return select_cortical_channels(experiment, subject, lfps)


def select_white_matter_channels(
    experiment: str,
    subject: str,
    da: xr.DataArray,
    wm_thickness: float = 200,
) -> xr.DataArray:
    s3 = get_project("shared")
    params = s3.load_experiment_subject_params(experiment, subject)
    probe = params["spwrProbe"]
    [hc_lo, hc_hi] = params["probes"][probe]["structureBounds"]["hippocampus"]
    is_white_matter = (da["y"] > hc_hi) & (da["y"] < (hc_hi + wm_thickness))
    return da.sel({"channel": is_white_matter})


def open_white_matter_lfps(subject: str, experiment: str, **kwargs) -> xr.DataArray:
    lfps = open_lfps(subject, experiment, **kwargs)
    return select_white_matter_channels(experiment, subject, lfps)


def get_hippocampal_subregion_channels(
    subject: str, experiment: str, da: xr.DataArray, key: str, radius: float = 80
) -> Tuple[np.ndarray, np.ndarray]:
    s3 = get_project("shared")
    subregions = s3.load_experiment_subject_json(
        experiment, subject, WNEFiles.HIPPOCAMPAL_SUBREGIONS
    )
    lo = subregions[key] - radius
    hi = subregions[key] + radius
    y = da["y"].load()
    mask = (y >= lo) & (y <= hi)
    ids = y["channel"][mask].values
    return mask, ids


# TODO: This doesn't need to be dedicated function.
def get_cortical_bandpower_file(
    subject: str, experiment: str, band: Bands
) -> pathlib.Path:
    return get_project("seahorse").get_experiment_subject_file(
        experiment, subject, Files.CORTICAL_BANDPOWER(band)
    )


# TODO: This doesn't need to be dedicated function.
def get_cortical_psds_file(subject: str, experiment: str) -> pathlib.Path:
    return get_project("seahorse").get_experiment_subject_file(
        experiment, subject, Files.CX_PSDS
    )


# TODO: This doesn't need to be dedicated function.
def get_hippocampal_bandpower_file(
    subject: str, experiment: str, band: Bands
) -> pathlib.Path:
    return get_project("seahorse").get_experiment_subject_file(
        experiment, subject, Files.HIPPOCAMPAL_BANDPOWER(band)
    )


# TODO: This doesn't need to be dedicated function.
def get_hippocampal_psds_file(subject: str, experiment: str) -> pathlib.Path:
    return get_project("seahorse").get_experiment_subject_file(
        experiment, subject, Files.HIPPOCAMPAL_PSDS
    )


# TODO: These is obsolete, and can be replaced with something much simpler.
# See ecephys.wne.sglx.utils.load_consolidated_artifacts().
def load_artifacts(
    experiment: str, subject: str, as_hypnogram: bool = True
) -> pd.DataFrame:
    """Stored, and returned, as float times."""
    s3 = get_project("shared")
    probe = get_spw_probe(experiment, subject)
    artifacts_file = s3.get_experiment_subject_file(
        experiment, subject, f"{probe}.lf.{WNEFiles.ARTIFACTS}"
    )
    if artifacts_file.is_file():
        artifacts = utils.read_htsv(artifacts_file)
        if as_hypnogram:
            if probe != "imec0":
                warnings.warn(
                    f"Attempting to load artifacts from {probe}. Is the hypnogram in the same timebase?"
                )
            artifacts["state"] = "Art"
    else:
        artifacts = pd.DataFrame(
            columns=["start_time", "end_time", "duration", "state"]
        )

    if as_hypnogram:
        artifacts = hypnogram.FloatHypnogram(artifacts)

    return artifacts


# TODO: These is obsolete, and can be replaced with something much simpler.
# It also probably doesn't need to be its own function.
def drop_artifacts(
    da: xr.DataArray, artifacts: hypnogram.FloatHypnogram, hg: hypnogram.FloatHypnogram
) -> xr.DataArray:
    is_art = artifacts.covers_time(da.time.values)
    is_art = is_art | hg.keep_states(ARTIFACT_STATES).covers_time(da.time.values)
    return da.drop_isel(time=is_art)


# TODO: Perhaps this should be in a ripples or localization module.
def get_relative_ripple_power_profile(lfp: xr.DataArray) -> xr.DataArray:
    lfp_spgs = xrsig.stft_psd(lfp)
    lfp_total_power = lfp_spgs.sum(dim="frequency").mean(dim="time")
    lfp_ripple_power = (
        lfp_spgs.sel(frequency=slice(*Bands.RIPPLE.value))
        .sum(dim="frequency")
        .mean(dim="time")
    )
    return lfp_ripple_power / lfp_total_power


# TODO: Perhaps belongs in an xrutils module.
def get_peaks_and_troughs(
    da: xr.DataArray, prominimence_pct=0.05, return_as: str = "y"
):
    assert len(da.shape) == 1, "Can only handle 1D data."

    def thresh(x):
        return float(np.abs(np.max(x) - np.min(x)) * prominimence_pct)

    peaks = scipy.signal.find_peaks(da, prominence=thresh(da))[0]
    troughs = scipy.signal.find_peaks(-da, prominence=thresh(-da))[0]
    return da[return_as].values[peaks], da[return_as].values[troughs]


def read_shifts(subject: str, experiment: str) -> pd.DataFrame:
    nb = get_project("seahorse")
    shift_info_file = nb.get_experiment_subject_file(
        experiment, subject, Files.RIPPLE_POWER_SHIFT
    )
    return (
        pd.read_parquet(shift_info_file)
        .sort_values("block_start_time")
        .reset_index(drop=True)
    )


def get_block_shifts(da: xr.DataArray, shifts: pd.DataFrame) -> np.ndarray:
    block_shifts = []
    for block in xrsig.iterate_timeseries_chunks(da):
        shifts_in_block = (shifts["block_start_time"] >= block["time"].values[0]) & (
            shifts["block_end_time"] <= block["time"].values[-1]
        )
        shifts_contain_block = (
            shifts["block_start_time"] <= block["time"].values[0]
        ) & (shifts["block_end_time"] >= block["time"].values[-1])
        block_shifts.append(
            int(
                shifts.loc[shifts_in_block | shifts_contain_block, "block_shift"].mean()
            )
        )
    return np.array(block_shifts).reshape(da.data.blocks.shape)


def get_shifted_channel_masks(
    shifts: pd.DataFrame, da: xr.DataArray, unshifted_channel_ids: np.ndarray
) -> xr.DataArray:
    unshifted_channel_mask = xr.where(
        da["channel"].isin(unshifted_channel_ids), True, False
    )
    return xr.concat(
        [
            unshifted_channel_mask.shift({"channel": shift}, fill_value=False)
            for shift in get_block_shifts(da, shifts).squeeze()
        ],
        dim="block",
    )


def get_event_coupling(
    evtsA: pd.DataFrame,
    evtsB: pd.DataFrame,
    prefixB: str,
    colsB: list[str] = ["zlog_amp", "zlog_duration"],
) -> pd.DataFrame:
    """For each event of type A, get the number of events of type B that overlap with it, and mean features of those events.
    Events must have start_time and end_time columns. Columns {prefixB}_count, and {prefixB}_{col} for col in `colsB` will be added to evtsA.
    """
    evtsA = evtsA.copy()
    evtsA[f"{prefixB}_count"] = 0
    for col in colsB:
        evtsA[f"{prefixB}_{col}"] = np.nan
    for ix in evtsA.index:
        a = evtsA.loc[ix]
        matches = (evtsB["end_time"] >= a["start_time"]) & (
            evtsB["start_time"] <= a["end_time"]
        )
        evtsA.loc[ix, f"{prefixB}_count"] = matches.sum()
        for col in colsB:
            evtsA.loc[ix, f"{prefixB}_{col}"] = evtsB.loc[matches, col].mean()
    return evtsA


def get_trial_shifts(da: xr.DataArray, shifts: pd.DataFrame) -> np.ndarray:
    trial_shifts = []
    for trial_time in da["event"].values:
        shifts_contain_trial = (shifts["block_start_time"] <= trial_time) & (
            shifts["block_end_time"] >= trial_time
        )
        trial_shifts.append(int(shifts.loc[shifts_contain_trial, "block_shift"].mean()))
    return np.array(trial_shifts)


def get_trial_shifted_channel_masks(
    shifts: pd.DataFrame, da: xr.DataArray, unshifted_channel_ids: np.ndarray
) -> xr.DataArray:
    unshifted_channel_mask = xr.where(
        da["channel"].isin(unshifted_channel_ids), True, False
    )
    return xr.concat(
        [
            unshifted_channel_mask.shift({"channel": shift}, fill_value=False)
            for shift in get_trial_shifts(da, shifts).squeeze()
        ],
        dim="event",
    ).assign_coords(event=da.event)
