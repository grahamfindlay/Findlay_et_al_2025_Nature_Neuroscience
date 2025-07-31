import warnings

import brainglobe_atlasapi
import numpy as np
import pandas as pd
import scipy.stats
import xarray as xr

import wisc_ecephys_tools as wet
from findlay2025a import agg, constants, core
from wisc_ecephys_tools import rats
from wisc_ecephys_tools.rats.constants import SleepDeprivationExperiments

_atlas = brainglobe_atlasapi.BrainGlobeAtlas("whs_sd_rat_39um")
_acronyms = ["Cx"] + _atlas.get_structure_descendants("Cx")


def agg_probe(experiment: str, subject: str, probe: str, iband: str) -> pd.DataFrame:
    nbsh = wet.get_sglx_project("seahorse")
    savefile = nbsh.get_experiment_subject_file(
        experiment, subject, f"{probe}.cx_{iband}.zarr"
    )
    da = xr.load_dataarray(savefile).dropna(dim="time")
    ds = xr.Dataset({iband: da})

    # TODO: Are these non-log-transformed z-scores ever used? Remove?
    ds[f"z{iband}"] = xr.apply_ufunc(scipy.stats.zscore, ds[iband])
    ds[f"log_{iband}"] = np.log10(ds[iband])
    ds[f"zlog_{iband}"] = xr.apply_ufunc(scipy.stats.zscore, ds[f"log_{iband}"])

    for acronym in ds["acronym"].values:
        if (ds[f"log_{iband}"].sel({"acronym": acronym}) < 0).any():
            warnings.warn(
                f"Low log_{iband} for {experiment} {subject} {probe} {acronym}. "
                "Please investigate."
            )

    hgs, _ = agg.get_hypnograms(experiment, subject)
    dfs = dict()
    for c, hg in hgs.items():
        c_mask = hg.covers_time(ds["time"])
        c_df = ds.isel({"time": c_mask}).to_dataframe().droplevel("time").reset_index()
        c_df["acronym"] = pd.Categorical(c_df["acronym"], categories=_acronyms)
        c_df["condition"] = pd.Categorical(
            [c] * len(c_df), categories=constants.CORE_STATISTICAL_CONDITIONS
        )
        dfs[c] = c_df
    return pd.concat(dfs.values())


def agg_project(iband: str) -> pd.DataFrame:
    available = rats.utils.get_subject_experiment_probe_tuples(
        experiment_filter=lambda x: x in SleepDeprivationExperiments
    )
    inproject = [
        (subj, str(expt)) for subj, expt in core.yield_subject_name_experiment_pairs()
    ]
    todo = [t for t in available if t[:2] in inproject]

    subjs = list(dict.fromkeys(t[0] for t in todo))
    expts = list(dict.fromkeys(t[1] for t in todo))
    prbs = list(dict.fromkeys(t[2] for t in todo))

    dfs = dict()
    for subject, experiment, probe in todo:
        print(f"Doing {subject} {experiment} {probe} {iband}")
        df = agg_probe(experiment, subject, probe, iband)
        df["subject"] = pd.Categorical([subject] * len(df), categories=subjs)
        df["experiment"] = pd.Categorical([experiment] * len(df), categories=expts)
        df["probe"] = pd.Categorical([probe] * len(df), categories=prbs)
        dfs[(subject, experiment, probe)] = df
    return pd.concat(dfs.values())


def do_project(iband: str):
    c_pwr = agg_project(iband)
    c_means = (
        c_pwr.groupby(
            ["subject", "experiment", "probe", "acronym", "condition"], observed=True
        )
        .mean(numeric_only=True)
        .add_prefix("cx_mean_")
    )  # Because observed=True, this will drop contrasts where a subject is missing either of the contrast conditions
    contrasts = agg.get_contrasts(c_means)

    nbsh = wet.get_sglx_project("seahorse")
    c_means.to_parquet(nbsh.get_project_file(f"cx_{iband}_condition_means.pqt"))
    contrasts.to_parquet(nbsh.get_project_file(f"cx_{iband}_condition_contrasts.pqt"))
