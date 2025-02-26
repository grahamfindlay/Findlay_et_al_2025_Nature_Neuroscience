import pandas as pd
import wisc_ecephys_tools as wet
from ecephys import hypnogram, wne

from findlay2025a import constants, core
from findlay2025a.constants import Experiments as Exps


def get_conveyor_over_water_period(
    experiment: str, wne_subject: wne.sglx.SGLXSubject
) -> tuple[float, float]:
    s3 = core.get_project("shared")
    params = s3.load_experiment_subject_params(experiment, wne_subject.name)
    probe = params["hypnogram_probe"]

    start = pd.to_datetime(params["conveyor_over_water_start"])
    start = wne_subject.dt2t(experiment, probe, start)

    end = pd.to_datetime(params["conveyor_over_water_end"])
    end = wne_subject.dt2t(experiment, probe, end)

    return (start, end)


def get_sleep_deprivation_period(
    experiment: str, wne_subject: wne.sglx.SGLXSubject
) -> tuple[float, float]:
    if experiment == Exps.NOD:
        return wet.shared.get_novel_objects_period(experiment, wne_subject)
    if experiment == Exps.COW:
        return get_conveyor_over_water_period(experiment, wne_subject)
    if experiment == Exps.CTN:
        start = min(get_conveyor_over_water_period(experiment, wne_subject))
        end = max(wet.shared.get_novel_objects_period(experiment, wne_subject))
        return start, end


def get_extended_wake_hypnogram(
    full_hg: hypnogram.FloatHypnogram,
    experiment: str,
    wne_subject: wne.sglx.SGLXSubject,
) -> hypnogram.FloatHypnogram:
    sd_start, sd_end = get_sleep_deprivation_period(experiment, wne_subject)
    five_minutes = pd.to_timedelta("5m").total_seconds()
    is_nod = (full_hg["start_time"] >= (sd_start - five_minutes)) & (
        full_hg["end_time"] <= sd_end
    )  # Note that this is a very strict criterion
    is_nodata = full_hg["state"] == "NoData"

    # Temporarily relabel NoData during SD, since we KNOW it is actually wake, and we want the algorithm that finds consolidated wake to account for this.
    full_hg.loc[is_nod & is_nodata, "state"] = "NoDataWake"
    matches = full_hg.get_consolidated(
        ["Wake", "Artifact", "Other", "NoDataWake"],
        minimum_time=(sd_end - sd_start) * 0.8,
        minimum_endpoint_bout_duration=pd.to_timedelta("120s").total_seconds(),
        maximum_antistate_bout_duration=pd.to_timedelta("90s").total_seconds(),
        frac=0.95,
    )
    full_hg.loc[is_nod & is_nodata, "state"] = "NoData"  # Restore NoData
    return matches[0].keep_states(["Wake"])


def get_conveyor_over_water_hypnogram(
    full_hg: hypnogram.FloatHypnogram,
    experiment: str,
    sglx_subject: wne.sglx.SGLXSubject,
) -> hypnogram.FloatHypnogram:
    (cow_start, cow_end) = get_conveyor_over_water_period(experiment, sglx_subject)
    return full_hg.trim(cow_start, cow_end)


def get_sleep_deprivation_wake_hypnogram(
    full_hg: hypnogram.FloatHypnogram,
    experiment: str,
    wne_subject: wne.sglx.SGLXSubject,
) -> hypnogram.FloatHypnogram:
    sd_start, sd_end = get_sleep_deprivation_period(experiment, wne_subject)
    return full_hg.trim(sd_start, sd_end).keep_states(["Wake"])


def get_post_deprivation_day2_light_period_hypnogram(
    full_hg: hypnogram.Hypnogram,
    experiment: str,
    sglx_subject: wne.sglx.SGLXSubject,
    sleep_deprivation_end: float,
) -> hypnogram.FloatHypnogram:
    d2lp_hg = wet.shared.get_day2_light_period_hypnogram(
        full_hg, experiment, sglx_subject
    )
    return d2lp_hg.trim(sleep_deprivation_end, d2lp_hg["end_time"].max())


def get_circadian_match_hypnogram(
    full_hg: hypnogram.FloatHypnogram, start: float, end: float
) -> hypnogram.FloatHypnogram:
    match_start = start - pd.to_timedelta("24h").total_seconds()
    match_end = end - pd.to_timedelta("24h").total_seconds()
    return full_hg.trim(match_start, match_end)


def compute_statistical_condition_hypnograms(
    lib_hg: hypnogram.FloatHypnogram,  # Full experiment hypnogram. "Liberal" because it provides the most accurate representation of the animal's sleep-wake state.
    cons_hg: hypnogram.FloatHypnogram,  # Full experiment hypnogram. "Conservative" because additional artifactual periods have been marked.
    experiment: str,
    sglx_subject: wne.sglx.SGLXSubject,
) -> dict[str, hypnogram.FloatHypnogram]:
    nrem_duration: float = pd.to_timedelta("1:00:00").total_seconds()
    wake_duration: float = pd.to_timedelta("1:00:00").total_seconds()
    rem_duration: float = pd.to_timedelta("0:10:00").total_seconds()

    hgs = dict()
    hgs["full_liberal"] = lib_hg
    hgs["full_conservative"] = cons_hg

    d1_hg = wet.shared.get_day1_hypnogram(cons_hg, experiment, sglx_subject)
    hgs["bsl_wake"] = d1_hg.keep_states(["Wake"])
    hgs["bsl_rem"] = d1_hg.keep_states(["REM"])

    d1lp_hg = wet.shared.get_day1_light_period_hypnogram(
        cons_hg, experiment, sglx_subject
    )
    hgs["early_bsl_nrem"] = d1lp_hg.keep_states(["NREM"]).keep_first(nrem_duration)
    hgs["early_bsl_rem"] = d1lp_hg.keep_states(["REM"]).keep_first(rem_duration)

    d1dp_hg = wet.shared.get_day1_dark_period_hypnogram(
        cons_hg, experiment, sglx_subject
    )
    hgs["last_bsl_nrem"] = d1dp_hg.keep_states(["NREM"]).keep_last(nrem_duration)
    hgs["last_bsl_rem"] = d1dp_hg.keep_states(["REM"]).keep_last(rem_duration)

    ewk_hg = get_extended_wake_hypnogram(lib_hg, experiment, sglx_subject)
    ewk_hg = cons_hg.trim(
        ewk_hg["start_time"].min(), ewk_hg["end_time"].max()
    ).keep_states(["Wake"])
    hgs["ext_wake"] = ewk_hg
    hgs["early_ext_wake"] = ewk_hg.keep_first(wake_duration)
    hgs["late_ext_wake"] = ewk_hg.keep_last(wake_duration)

    sd_hg = get_sleep_deprivation_wake_hypnogram(lib_hg, experiment, sglx_subject)
    sd_hg = cons_hg.trim(
        sd_hg["start_time"].min(), sd_hg["end_time"].max()
    ).keep_states(["Wake"])
    hgs["sd_wake"] = sd_hg
    hgs["early_sd_wake"] = sd_hg.keep_first(wake_duration)
    hgs["late_sd_wake"] = sd_hg.keep_last(wake_duration)

    if experiment in [Exps.NOD, Exps.CTN]:
        nod_hg = wet.shared.get_novel_objects_hypnogram(
            cons_hg, experiment, sglx_subject
        ).keep_states(["Wake"])
        hgs["nod_wake"] = nod_hg
        hgs["early_nod_wake"] = nod_hg.keep_first(wake_duration)
        hgs["late_nod_wake"] = nod_hg.keep_last(wake_duration)

    if experiment in [Exps.COW, Exps.CTN]:
        cow_hg = get_conveyor_over_water_hypnogram(
            cons_hg, experiment, sglx_subject
        ).keep_states(["Wake"])
        hgs["cow_wake"] = cow_hg
        hgs["early_cow_wake"] = cow_hg.keep_first(wake_duration)
        hgs["late_cow_wake"] = cow_hg.keep_last(wake_duration)

    if experiment == Exps.CTN:
        hgs["ctn_wake"] = hgs["sd_wake"]
        hgs["early_ctn_wake"] = hgs["early_sd_wake"]
        hgs["late_ctn_wake"] = hgs["late_sd_wake"]

    pdd2lp_hg = get_post_deprivation_day2_light_period_hypnogram(
        cons_hg,
        experiment,
        sglx_subject,
        sleep_deprivation_end=ewk_hg["end_time"].max(),
    )
    hgs["early_rec_nrem"] = pdd2lp_hg.keep_states(["NREM"]).keep_first(nrem_duration)
    hgs["early_rec_rem"] = pdd2lp_hg.keep_states(["REM"]).keep_first(rem_duration)
    hgs["late_rec_nrem"] = pdd2lp_hg.keep_states(["NREM"]).keep_last(nrem_duration)
    hgs["late_rec_rem"] = pdd2lp_hg.keep_states(["REM"]).keep_last(rem_duration)

    thirty_minutes = pd.to_timedelta("30m").total_seconds()
    hgs["early_rec_nrem_match"] = get_circadian_match_hypnogram(
        cons_hg,
        start=hgs["early_rec_nrem"]["start_time"].min() - thirty_minutes,
        end=hgs["early_rec_nrem"]["end_time"].max() + thirty_minutes,
    ).keep_states(["NREM"])
    hgs["early_rec_rem_match"] = get_circadian_match_hypnogram(
        cons_hg,
        start=hgs["early_rec_rem"]["start_time"].min() - thirty_minutes,
        end=hgs["early_rec_rem"]["end_time"].max() + thirty_minutes,
    ).keep_states(["REM"])

    d2dp_hg = wet.shared.get_day2_dark_period_hypnogram(
        cons_hg, experiment, sglx_subject
    )
    hgs["last_rec_nrem"] = d2dp_hg.keep_states(["NREM"]).keep_last(nrem_duration)
    hgs["last_rec_rem"] = d2dp_hg.keep_states(["REM"]).keep_last(rem_duration)

    return hgs


def load_statistical_condition_hypnograms(
    experiment: str, subject: str, include_full_conservative: bool = False
) -> dict[str, hypnogram.FloatHypnogram]:
    nb = core.get_project("seahorse")
    hgs = {}
    conditions_to_load = (
        constants.ALL_STATISTICAL_CONDITIONS + ["full_conservative"]
        if include_full_conservative
        else constants.ALL_STATISTICAL_CONDITIONS
    )
    for c in conditions_to_load:
        try:
            hgs[c] = hypnogram.FloatHypnogram.from_htsv(
                nb.get_experiment_subject_file(
                    experiment, subject, f"{c}_hypnogram.htsv"
                )
            )
        except FileNotFoundError:
            pass
    return hgs
