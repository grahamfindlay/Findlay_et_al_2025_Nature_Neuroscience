from collections import defaultdict
from functools import lru_cache

import brainglobe_atlasapi
import pandas as pd

import wisc_ecephys_tools as wet
from ecephys import units, wne
from findlay2025a import core


def get_threshold_kwargs():
    return dict(
        all=dict(
            required_threshold="all",
            isolation_threshold=None,
            false_negatives_threshold=None,
            presence_threshold=None,
        ),
        mua=dict(
            required_threshold="conservative",
            isolation_threshold=None,
            false_negatives_threshold=None,
            presence_threshold=None,
        ),
        sua_permissive=dict(
            required_threshold="conservative",
            isolation_threshold="permissive",
            false_negatives_threshold="permissive",
            presence_threshold=None,
        ),
        sua_moderate=dict(
            required_threshold="conservative",
            isolation_threshold="moderate",
            false_negatives_threshold="moderate",
            presence_threshold=None,
        ),
        sua_conservative=dict(
            required_threshold="conservative",
            isolation_threshold="conservative",
            false_negatives_threshold="conservative",
            presence_threshold=None,
        ),
    )


@lru_cache(maxsize=4)
def get_nod_sortings(
    spw_probe_only: bool = False, expand_probes: bool = False
) -> list[tuple[str, ...]]:
    nod = wet.rats.constants.SleepDeprivationExperiments.NOD
    s3 = wet.projects.get_wne_project("shared")
    sortings = wet.rats.utils.get_subject_experiment_probe_tuples(
        experiment_filter=lambda x: x == nod
    )
    sortings = [
        (s, p)
        for s, e, p in sortings
        if wet.rats.utils.has_sorting(s, e, p, s3)
        and wet.rats.utils.has_anatomy(s, e, p, s3)
        and wet.rats.utils.has_hypnogram(s, e, None, s3)
    ]
    spw_subjects = [
        sub for sub, _ in core.yield_subject_name_experiment_pairs(experiments=[nod])
    ]
    sortings = [(s, p) for (s, p) in sortings if s in spw_subjects]
    if spw_probe_only:
        sortings = [(s, p) for (s, p) in sortings if core.get_spw_probe(nod, s) == p]
    if not expand_probes:
        s2p = defaultdict(list)
        for s, p in sortings:
            s2p[s].append(p)
        sortings = [(s, tuple(p)) for s, p in s2p.items()]
    return sortings


def load_nod_multiprobe_sorting(subject: str, **threshold_kwargs) -> units.MultiSIKS:
    s3 = core.get_project("shared")
    experiment = wet.rats.constants.SleepDeprivationExperiments.NOD

    units_ready = {subj: prbs for (subj, prbs) in get_nod_sortings()}
    probes = units_ready[subject]
    mps = wne.sglx.siutils.load_multiprobe_sorting(
        s3,
        subject,
        experiment,
        probes=probes,
        wneAnatomyProject=s3,
    )

    simple_filters, callable_filters = wne.siutils.get_quality_metric_filters(
        **threshold_kwargs
    )
    simple_filters_by_probe = {probe: simple_filters for probe in probes}
    callable_filters_by_probe = {probe: callable_filters for probe in probes}
    mps = mps.refine_clusters(
        simple_filters_by_probe, callable_filters_by_probe, include_nans=True
    )

    return mps


def load_nod_spw_probe_sorting(
    subject: str, **threshold_kwargs
) -> units.SpikeInterfaceKilosortSorting:
    s3 = core.get_project("shared")
    experiment = wet.rats.constants.SleepDeprivationExperiments.NOD

    spw_probe = core.get_spw_probe(experiment, subject)
    s = wne.sglx.siutils.load_singleprobe_sorting(
        s3,
        subject,
        experiment,
        spw_probe,
        wneAnatomyProject=s3,
        allow_no_sync_file=False,
    )

    simple_filters, callable_filters = wne.siutils.get_quality_metric_filters(
        **threshold_kwargs
    )
    s = s.refine_clusters(simple_filters, callable_filters, include_nans=True)

    return s


def hippocampus_to_waxholm(acronym: str) -> str:
    if "CA1" in acronym:
        return "CA1"
    elif "CA2" in acronym:
        return "CA2"
    elif "CA3" in acronym:
        return "CA3"
    elif "DG" in acronym:
        return "DG"
    else:
        return acronym


def add_major_regions(df: pd.DataFrame, as_booleans: bool = True) -> pd.DataFrame:
    """Use the `acronym` column of a dataframe to add major region information.

    Assumes the `acronym` column contains Waxholm atlas regions and adds major
    region columns (if `as_booleans == True`) or labels to a new `region` column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with an `acronym` column containing Waxholm atlas region names.
    as_booleans : bool, default True
        If True, adds boolean columns for each major region (is_hippocampal,
        is_thalamic, is_cortical, is_other). If False, adds a single `region`
        column with the region type (hippocampal, thalamic, cortical, other).

    Returns
    -------
    pd.DataFrame
        DataFrame with added region information.
    """
    atlas = brainglobe_atlasapi.BrainGlobeAtlas("whs_sd_rat_39um")
    df["is_hippocampal"] = df["acronym"].isin(
        atlas.get_structure_descendants("HF") + ["HF"]
    )
    df["is_thalamic"] = df["acronym"].isin(
        atlas.get_structure_descendants("Thal-D") + ["Thal-D"]
    )
    df["is_cortical"] = ~df["is_hippocampal"] & df["acronym"].isin(
        atlas.get_structure_descendants("Cx") + ["Cx"]
    )
    df["is_other"] = ~df["is_hippocampal"] & ~df["is_thalamic"] & ~df["is_cortical"]

    if not as_booleans:
        # Make sure the region columns are mutually exclusive
        region_cols = ["is_hippocampal", "is_thalamic", "is_cortical", "is_other"]
        assert (df[region_cols].sum(axis=1) == 1).all()

        # Melt the region columns into a single 'region' column
        df_melted = df.melt(
            id_vars=[col for col in df.columns if col not in region_cols],
            value_vars=region_cols,
            var_name="region_flag",
            value_name="is_region",
        )

        # Keep only rows where is_region is True and clean up the region names
        df = df_melted[df_melted["is_region"]].copy()
        df["region"] = df["region_flag"].str.replace("is_", "")
        df = df.drop(columns=["region_flag", "is_region"])

    return df
