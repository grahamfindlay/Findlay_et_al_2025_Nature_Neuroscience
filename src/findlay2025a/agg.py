import itertools
import warnings
from types import MappingProxyType
from typing import Callable, Tuple

import numpy as np
import pandas as pd
import scipy.stats
import xarray as xr

from ecephys.hypnogram import FloatHypnogram
from ecephys.wne.constants import Files as wneFiles
from findlay2025a import core, hypnograms
from findlay2025a.constants import (
    CORE_STATISTICAL_CONDITIONS,
    COW_STATISTICAL_CONDITIONS,
    CTN_STATISTICAL_CONDITIONS,
    NOD_STATISTICAL_CONDITIONS,
    Bands,
    Experiments,
)

CONTRASTS = MappingProxyType(
    {
        "nrem_rebound": ("early_rec_nrem", "early_rec_nrem_match"),
        "nrem_surge": ("early_rec_nrem", "early_bsl_nrem"),
        "nrem_rec_decline": ("early_rec_nrem", "late_rec_nrem"),
        "nrem_bsl_decline": ("early_bsl_nrem", "early_rec_nrem_match"),
        "ext_wake_incline": ("late_ext_wake", "early_ext_wake"),
        "sd_wake_incline": ("early_sd_wake", "late_sd_wake"),
    }
)


def aggregated_events_wide_to_long(
    evts_wide: pd.DataFrame, condition_cols: list[str] = None, drop: bool = True
) -> pd.DataFrame:
    if condition_cols is None:
        condition_cols = set(CORE_STATISTICAL_CONDITIONS)
    evts_long = pd.DataFrame()
    for c in condition_cols:
        c_evts = evts_wide[evts_wide[c]].copy()
        c_evts["condition"] = c
        evts_long = pd.concat([evts_long, c_evts], ignore_index=True)
    if drop:
        evts_long = evts_long.drop(columns=condition_cols)
    return evts_long


def get_contrasts(
    c_means: pd.DataFrame,
    c_sums: pd.DataFrame = None,
    c_rates: pd.DataFrame = None,
    rate_columns: list[str] = None,
) -> pd.DataFrame:
    contrast_inputs = c_means
    if c_sums is not None:
        contrast_inputs = contrast_inputs.join(c_sums)
    if c_rates is not None:
        contrast_inputs = contrast_inputs.join(c_rates[rate_columns])

    contrast_dfs = {}
    for contrast, (condition_a, condition_b) in CONTRASTS.items():
        df = contrast_inputs.xs(condition_a, level="condition") - contrast_inputs.xs(
            condition_b, level="condition"
        )
        contrast_dfs[contrast] = df.add_suffix(f"_{contrast}")
    contrasts = pd.concat(list(contrast_dfs.values()), axis="columns")
    return contrasts


def aggregate_cortical_bandpowers() -> pd.DataFrame:
    def _load_cortical_bandpower(
        subject: str, experiment: str, band: Bands
    ) -> xr.DataArray:
        return xr.load_dataarray(
            core.get_cortical_bandpower_file(subject, experiment, band)
        )

    return _aggregate_bandpowers(_load_cortical_bandpower)


def aggregate_hippocampal_bandpowers() -> pd.DataFrame:
    def _load_hippocampal_bandpower(
        subject: str, experiment: str, band: Bands
    ) -> xr.DataArray:
        return xr.load_dataarray(
            core.get_hippocampal_bandpower_file(subject, experiment, band)
        )

    return _aggregate_bandpowers(_load_hippocampal_bandpower)


def get_hypnograms(experiment: str, subject: str) -> FloatHypnogram:
    hgs = hypnograms.load_statistical_condition_hypnograms(
        experiment, subject, include_full_conservative=True
    )
    omit = core.ARTIFACT_STATES + ["NoData"]
    fc = hgs.pop("full_conservative").drop_states(omit)
    return hgs, fc


def _aggregate_bandpowers(
    loader_function: Callable[[str, str, str], xr.DataArray],
) -> pd.DataFrame:
    res = pd.DataFrame()
    for subject, experiment in core.yield_subject_name_experiment_pairs():
        hgs, hg = get_hypnograms(experiment, subject)

        delta = loader_function(subject, experiment, Bands.DELTA)
        theta = loader_function(subject, experiment, Bands.THETA)
        gamma = loader_function(subject, experiment, Bands.GAMMA)
        eta = loader_function(subject, experiment, Bands.VLAD)
        ds = xr.Dataset({"delta": delta, "theta": theta, "gamma": gamma, "eta": eta})
        ds = ds.isel(time=hg.covers_time(ds["time"]))

        ds["tdr"] = ds["theta"] / ds["delta"]
        # TODO: Are these non-log-transformed z-scores ever used? Remove?
        ds["zdelta"] = xr.apply_ufunc(scipy.stats.zscore, ds["delta"])
        ds["ztheta"] = xr.apply_ufunc(scipy.stats.zscore, ds["theta"])
        ds["zgamma"] = xr.apply_ufunc(scipy.stats.zscore, ds["gamma"])
        ds["zeta"] = xr.apply_ufunc(scipy.stats.zscore, ds["eta"])
        ds["log_delta"] = np.log10(ds["delta"])
        ds["log_theta"] = np.log10(ds["theta"])
        ds["log_gamma"] = np.log10(ds["gamma"])
        ds["log_eta"] = np.log10(ds["eta"])
        ds["zlog_delta"] = xr.apply_ufunc(scipy.stats.zscore, ds["log_delta"])
        ds["zlog_theta"] = xr.apply_ufunc(scipy.stats.zscore, ds["log_theta"])
        ds["zlog_gamma"] = xr.apply_ufunc(scipy.stats.zscore, ds["log_gamma"])
        ds["zlog_eta"] = xr.apply_ufunc(scipy.stats.zscore, ds["log_eta"])

        for var in ["log_delta", "log_theta", "log_gamma", "log_eta"]:
            if (ds[var] < 0).any():
                warnings.warn(
                    f"Low {var} for {subject} in {experiment}. Please investigate.",
                    UserWarning,
                )

        df = ds.to_dataframe().reset_index()
        df["subject"] = subject
        df["experiment"] = experiment

        for c, hg in hgs.items():
            df[c] = hg.covers_time(df["time"])
        res = pd.concat([res, df], axis=0, ignore_index=True)
    for c in CORE_STATISTICAL_CONDITIONS:
        res[c] = res[c].fillna(False)
    return res.drop(columns=["time"])


def aggregate_emg() -> pd.DataFrame:
    omit = core.ARTIFACT_STATES + ["NoData"]
    res = pd.DataFrame()
    for subject, experiment in core.yield_subject_name_experiment_pairs():
        hgs = hypnograms.load_statistical_condition_hypnograms(
            experiment, subject, include_full_conservative=True
        )
        hg = hgs.pop("full_conservative").drop_states(omit)

        emg = xr.load_dataarray(
            core.get_project("seahorse").get_experiment_subject_file(
                experiment, subject, wneFiles.EMG
            )
        )
        emg = emg.drop_duplicates(dim="time", keep="first").sortby("time")
        da = xr.Dataset({"emg": emg})
        da = da.isel(time=hg.covers_time(da["time"]))

        da["zemg"] = xr.apply_ufunc(scipy.stats.zscore, da["emg"])
        da["log_emg"] = np.log10(da["emg"])
        da["zlog_emg"] = xr.apply_ufunc(scipy.stats.zscore, da["log_emg"])
        df = da.to_dataframe().reset_index()
        df["subject"] = subject
        df["experiment"] = experiment

        for c, hg in hgs.items():
            df[c] = hg.covers_time(df["time"])
        res = pd.concat([res, df], axis=0, ignore_index=True)
    for c in CORE_STATISTICAL_CONDITIONS:
        res[c] = res[c].fillna(False)
    return res.drop(columns=["time"])


def get_experiment_dataframe() -> pd.DataFrame:
    pairs = list(core.yield_subject_name_experiment_pairs())
    return pd.DataFrame(pairs, columns=["subject", "experiment"])


def get_state_dataframe(
    states: list[str] = ["Wake", "NREM", "IS", "REM"],
) -> pd.DataFrame:
    rows = [
        (*pair, state)
        for pair, state in itertools.product(
            core.yield_subject_name_experiment_pairs(), states
        )
    ]
    return pd.DataFrame(rows, columns=["subject", "experiment", "state"])


def get_condition_dataframe() -> pd.DataFrame:
    nod = [
        (*pair, condition)
        for pair, condition in itertools.product(
            core.yield_subject_name_experiment_pairs([Experiments.NOD]),
            NOD_STATISTICAL_CONDITIONS,
        )
    ]
    cow = [
        (*pair, condition)
        for pair, condition in itertools.product(
            core.yield_subject_name_experiment_pairs([Experiments.COW]),
            COW_STATISTICAL_CONDITIONS,
        )
    ]
    ctn = [
        (*pair, condition)
        for pair, condition in itertools.product(
            core.yield_subject_name_experiment_pairs([Experiments.CTN]),
            CTN_STATISTICAL_CONDITIONS,
        )
    ]
    return pd.DataFrame(nod + cow + ctn, columns=["subject", "experiment", "condition"])


def aggregate_events(
    loader_function: Callable[[str, str], pd.DataFrame],
    time_col: str,
    keep_cols: list[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    omit = core.ARTIFACT_STATES + ["NoData"]
    aggd_evts = pd.DataFrame()
    aggd_experiment_rates = (
        get_experiment_dataframe().set_index(["subject", "experiment"]).sort_index()
    )
    aggd_state_rates = (
        get_state_dataframe().set_index(["subject", "experiment", "state"]).sort_index()
    )
    aggd_condition_rates = (
        get_condition_dataframe()
        .set_index(["subject", "experiment", "condition"])
        .sort_index()
    )

    for subject, experiment in core.yield_subject_name_experiment_pairs():
        # Load the events, and label them with the subject and experiment.
        _evts = loader_function(subject, experiment)
        _evts["subject"] = subject
        _evts["experiment"] = experiment

        # Load the conditions and their hypnograms.
        hgs = hypnograms.load_statistical_condition_hypnograms(
            experiment, subject, include_full_conservative=True
        )
        # Drop states that we do not want to contribute to global rate estimation.
        hg = hgs.pop("full_conservative").drop_states(omit)
        _evts = _evts.iloc[hg.covers_time(_evts[time_col])]

        # Get event counts and durations. First for whole experiments...
        aggd_experiment_rates.loc[(subject, experiment), "count"] = len(_evts)
        aggd_experiment_rates.loc[(subject, experiment), "duration"] = hg[
            "duration"
        ].sum()

        # Then for states...
        states = list(aggd_state_rates.groupby("state").groups.keys())
        for state in states:
            hg_state = hg.keep_states([state])
            aggd_state_rates.loc[(subject, experiment, state), "count"] = (
                hg_state.covers_time(_evts[time_col]).sum()
            )
            aggd_state_rates.loc[(subject, experiment, state), "duration"] = hg_state[
                "duration"
            ].sum()

        # Label each event with the conditions it falls under.
        for condition, hg in hgs.items():
            _evts[condition] = hg.covers_time(_evts[time_col])
            aggd_condition_rates.loc[(subject, experiment, condition), "count"] = _evts[
                condition
            ].sum()
            aggd_condition_rates.loc[(subject, experiment, condition), "duration"] = hg[
                "duration"
            ].sum()
        aggd_evts = pd.concat([aggd_evts, _evts], axis=0, ignore_index=True)

    # Compute rates
    aggd_experiment_rates["rate"] = (
        aggd_experiment_rates["count"] / aggd_experiment_rates["duration"]
    )
    aggd_state_rates["rate"] = aggd_state_rates["count"] / aggd_state_rates["duration"]
    aggd_state_rates["rate_rel2total"] = (
        aggd_state_rates["rate"] / aggd_experiment_rates["rate"]
    )
    aggd_condition_rates["rate"] = (
        aggd_condition_rates["count"] / aggd_condition_rates["duration"]
    )
    aggd_condition_rates["rate_rel2total"] = (
        aggd_condition_rates["rate"] / aggd_experiment_rates["rate"]
    )

    # Convert any remaining NaNs to False.
    for condition in CORE_STATISTICAL_CONDITIONS:
        aggd_evts[condition] = aggd_evts[condition].fillna(False)

    # If the user only requested certain event properties, discard the others, to reduce data size
    if keep_cols is not None:
        keep_cols = ["subject", "experiment"] + keep_cols + CORE_STATISTICAL_CONDITIONS
    else:
        keep_cols = aggd_evts.columns

    return (
        aggd_evts[keep_cols],
        aggd_experiment_rates,
        aggd_state_rates,
        aggd_condition_rates,
    )
