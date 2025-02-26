import itertools
import warnings
from types import MappingProxyType
from typing import Callable, Tuple

import numpy as np
import pandas as pd
import scipy.stats
import xarray as xr
from ecephys.wne.constants import Files as wneFiles

from findlay2025a import core, hypnograms
from findlay2025a.constants import (
    ALL_STATISTICAL_CONDITIONS,
    COW_STATISTICAL_CONDITIONS,
    CTN_STATISTICAL_CONDITIONS,
    NOD_STATISTICAL_CONDITIONS,
    Bands,
    Experiments,
)

# TODO: Some of these are never used. Prune them.
CONTRASTS = MappingProxyType(
    {
        "nrem_rebound": ("early_rec_nrem", "early_rec_nrem_match"),
        "rem_rebound": ("early_rec_rem", "early_rec_rem_match"),
        "nrem_surge": ("early_rec_nrem", "early_bsl_nrem"),
        "rem_surge": ("early_rec_rem", "early_bsl_rem"),
        "nrem_rec_decline": ("early_rec_nrem", "late_rec_nrem"),
        "rem_rec_decline": ("early_rec_rem", "late_rec_rem"),
        "nrem_bsl_decline": ("early_bsl_nrem", "early_rec_nrem_match"),
        "rem_bsl_decline": ("early_bsl_rem", "early_rec_rem_match"),
        "ext_wake_incline": ("late_ext_wake", "early_ext_wake"),
        "sd_wake_incline": ("late_sd_wake", "early_sd_wake"),
        "nod_wake_incline": ("late_nod_wake", "early_nod_wake"),
        "cow_wake_incline": ("late_cow_wake", "early_cow_wake"),
    }
)


def aggregated_events_wide_to_long(
    evts_wide: pd.DataFrame, condition_cols: list[str] = None, drop: bool = True
) -> pd.DataFrame:
    if condition_cols is None:
        condition_cols = set(ALL_STATISTICAL_CONDITIONS)
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
        df = (
            contrast_inputs.loc[(slice(None), slice(None), condition_a)]
            - contrast_inputs.loc[(slice(None), slice(None), condition_b)]
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


def _aggregate_bandpowers(
    loader_function: Callable[[str, str, str], xr.DataArray],
) -> pd.DataFrame:
    omit = core.ARTIFACT_STATES + ["NoData"]
    res = pd.DataFrame()
    for subject, experiment in core.yield_subject_name_experiment_pairs():
        hgs = hypnograms.load_statistical_condition_hypnograms(
            experiment, subject, include_full_conservative=True
        )
        hg = hgs.pop("full_conservative").drop_states(omit)

        delta = loader_function(subject, experiment, Bands.DELTA)
        theta = loader_function(subject, experiment, Bands.THETA)
        gamma = loader_function(subject, experiment, Bands.GAMMA)
        da = xr.Dataset({"delta": delta, "theta": theta, "gamma": gamma})
        da = da.isel(time=hg.covers_time(da["time"]))

        da["tdr"] = da["theta"] / da["delta"]
        da["zdelta"] = scipy.stats.zscore(da["delta"])
        da["ztheta"] = scipy.stats.zscore(da["theta"])
        da["zgamma"] = scipy.stats.zscore(da["gamma"])
        da["log_delta"] = np.log10(da["delta"])
        da["log_theta"] = np.log10(da["theta"])
        da["log_gamma"] = np.log10(da["gamma"])
        da["zlog_delta"] = scipy.stats.zscore(da["log_delta"])
        da["zlog_theta"] = scipy.stats.zscore(da["log_theta"])
        da["zlog_gamma"] = scipy.stats.zscore(da["log_gamma"])

        for var in ["log_delta", "log_theta", "log_gamma"]:
            if (da[var] < 0).any():
                warnings.warn(
                    f"Low {var} for {subject} in {experiment}. Please investigate.",
                    UserWarning,
                )

        df = da.to_dataframe().reset_index()
        df["subject"] = subject
        df["experiment"] = experiment

        for c, hg in hgs.items():
            df[c] = hg.covers_time(df["time"])
        res = pd.concat([res, df], axis=0, ignore_index=True)
    for c in ALL_STATISTICAL_CONDITIONS:
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

        da["zemg"] = scipy.stats.zscore(da["emg"])
        da["log_emg"] = np.log10(da["emg"])
        da["zlog_emg"] = scipy.stats.zscore(da["log_emg"])
        df = da.to_dataframe().reset_index()
        df["subject"] = subject
        df["experiment"] = experiment

        for c, hg in hgs.items():
            df[c] = hg.covers_time(df["time"])
        res = pd.concat([res, df], axis=0, ignore_index=True)
    for c in ALL_STATISTICAL_CONDITIONS:
        res[c] = res[c].fillna(False)
    return res.drop(columns=["time"])


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
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    omit = core.ARTIFACT_STATES + ["NoData"]
    aggd_evts = pd.DataFrame()
    aggd_condition_rates = (
        get_condition_dataframe()
        .set_index(["subject", "experiment", "condition"])
        .sort_index()
    )
    aggd_experiment_rates = (
        get_condition_dataframe()[["subject", "experiment"]]
        .drop_duplicates(ignore_index=True)
        .set_index(["subject", "experiment"])
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

        # Get the total number of events outside the omitted states.
        aggd_experiment_rates.loc[(subject, experiment), "count"] = len(_evts)
        # Get the total amount of time spent outside the omitted states.
        aggd_experiment_rates.loc[(subject, experiment), "duration"] = hg[
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
    aggd_condition_rates["rate"] = (
        aggd_condition_rates["count"] / aggd_condition_rates["duration"]
    )
    aggd_experiment_rates["rate"] = (
        aggd_experiment_rates["count"] / aggd_experiment_rates["duration"]
    )
    aggd_condition_rates["rate_rel2total"] = (
        aggd_condition_rates["rate"] / aggd_experiment_rates["rate"]
    )

    # Convert any remaining NaNs to False.
    for condition in ALL_STATISTICAL_CONDITIONS:
        aggd_evts[condition] = aggd_evts[condition].fillna(False)

    # If the user only requested certain event properties, discard the others, to reduce data size
    if keep_cols is not None:
        keep_cols = ["subject", "experiment"] + keep_cols + ALL_STATISTICAL_CONDITIONS
    else:
        keep_cols = aggd_evts.columns

    return aggd_evts[keep_cols], aggd_experiment_rates, aggd_condition_rates
