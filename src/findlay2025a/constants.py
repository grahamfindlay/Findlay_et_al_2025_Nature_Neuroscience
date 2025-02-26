import math
from enum import Enum, StrEnum
from typing import Final, Tuple


class Experiments(StrEnum):
    NOD = "novel_objects_deprivation"
    COW = "conveyor_over_water"
    CTN = "conveyor_then_novelty"


class Bands(Enum):
    DELTA = (0.5, 4)
    THETA = (5, 10)
    SIGMA = (9, 16)
    GAMMA = (30, 150)
    SLOW_GAMMA = (30, 50)
    FAST_GAMMA = (50, 150)
    RIPPLE = (125, 250)
    MUA = (500, math.inf)

    AERY_JONES_SLOW_GAMMA = (20, 50)
    AERY_JONES_FAST_GAMMA = (50, 110)

    @property
    def low(self) -> float:
        return self.value[0]

    @property
    def high(self) -> float:
        return self.value[1]

    def contains(self, freq: float) -> bool:
        return self.low <= freq <= self.high


class Files(StrEnum):
    HIPPOCAMPAL_KCSD = "hippocampal_kcsd.zarr"
    KCSD_PARAMS = "kcsd_params.npz"

    SPW_PARAMS = "spw_params.npz"
    ESTM_SPWS = "estm_period_spws.pqt"
    SPW_DETECTION_CHANS = "classic_spw_detection_channels.nc"
    CLASSIC_SPWS = "classic_spws.pqt"  # TODO: Remove CLASSIC_ prefix
    POSTPROCESSED_SPWS = "postprocessed_spws.pqt"
    SPW_CONDITION_RATES = "spw_condition_rates.pqt"
    SPW_CONDITION_MEANS = "spw_condition_means.pqt"
    SPW_CONDITION_SUMS = "spw_condition_sums.pqt"
    SPW_CONDITION_CONTRASTS = "spw_condition_contrasts.pqt"
    SPW_LFP_SNIPPETS = "spw_lfp_snippets.nc"
    SPW_TRIGGERED_RIPPLE_POWER = "spw_triggered_ripple_power.nc"
    SPW_TRIGGERED_RIPPLE_FREQ = "spw_triggered_ripple_freq.nc"

    RIPPLE_PARAMS = "ripple_params.npz"
    ESTM_RIPPLES = "estm_period_ripples.pqt"
    RIPPLE_DETECTION_CHANS = "ripple_detection_channels.nc"
    RIPPLES = "ripples.pqt"
    POSTPROCESSED_RIPPLES = "postprocessed_ripples.pqt"
    RIPPLE_CONDITION_RATES = "ripple_condition_rates.pqt"
    RIPPLE_CONDITION_MEANS = "ripple_condition_means.pqt"
    RIPPLE_CONDITION_SUMS = "ripple_condition_sums.pqt"
    RIPPLE_CONDITION_CONTRASTS = "ripple_condition_contrasts.pqt"
    RIPPLE_LFP_SNIPPETS = "ripple_lfp_snippets.nc"
    RIPPLE_CSD_SNIPPETS = "ripple_csd_snippets.nc"
    RIPPLE_TRIGGERED_SPW_AMP = "ripple_triggered_spw_amp.nc"
    RIPPLE_FREQ = "ripple_freq.nc"
    RIPPLE_POWER = "ripple_power.nc"

    DSPK_PARAMS = "dspk_params.npz"
    ESTM_DSPKS = "estm_period_dspks.pqt"
    TMP_DPKS = "tmp_dpks.pqt"
    DPKS = "dpks.pqt"
    DSPKS = "dspks.pqt"
    POSTPROCESSED_DSPKS = "postprocessed_dspks.pqt"
    DSPK_CONDITION_RATES = "dspk_condition_rates.pqt"
    DSPK_CONDITION_MEANS = "dspk_condition_means.pqt"
    DSPK_CONDITION_SUMS = "dspk_condition_sums.pqt"
    DSPK_CONDITION_CONTRASTS = "dspk_condition_contrasts.pqt"

    ESTM_PERIOD_KCSD_PROFILE = "estm_period_kcsd_profile.nc"
    ESTM_PERIOD_RIPPLE_PROFILE = "estm_period_ripple_power_profile.nc"
    RIPPLE_POWER_SHIFT = "ripple_power_shift.pqt"

    CX_PSDS = "cx_psds_by_condition.nc"
    HIPPOCAMPAL_PSDS = "hipp_psds_by_condition.nc"

    CX_BANDPOWER_MEANS = "cx_bandpower_condition_means.pqt"
    CX_BANDPOWER_SUM = "cx_bandpower_condition_sums.pqt"
    CX_BANDPOWER_CONTRASTS = "cx_bandpower_condition_contrasts.pqt"

    HIPPOCAMPAL_BANDPOWER_MEANS = "hipp_bandpower_condition_means.pqt"
    HIPPOCAMPAL_BANDPOWER_SUM = "hipp_bandpower_condition_sums.pqt"
    HIPPOCAMPAL_BANDPOWER_CONTRASTS = "hipp_bandpower_condition_contrasts.pqt"

    EMG_CONDITION_MEANS = "emg_condition_means.pqt"
    EMG_CONDITION_SUMS = "emg_condition_sums.pqt"
    EMG_CONDITION_CONTRASTS = "emg_condition_contrasts.pqt"

    COMBINED_CONDITION_MEASURES = "combined_condition_measures.pqt"
    COMBINED_CONDITION_CONTRASTS = "combined_condition_contrasts.pqt"

    @staticmethod
    def HIPPOCAMPAL_BANDPOWER(band: Bands) -> str:
        return f"hipp_{band.name.lower()}.nc"

    @staticmethod
    def CORTICAL_BANDPOWER(band: Bands) -> str:
        return f"cx_{band.name.lower()}.nc"

    @staticmethod
    def COMODULOMGRAMS(state: str) -> str:
        return f"{state}_comodulograms.nc"

    @staticmethod
    def COMODULOMGRAM_SURROGATE_MAXES(state: str) -> str:
        return f"{state}_comodulogram_surrogate_maxes.nc"


# TODO: Most of these are never used. Prune them.
ALL_STATISTICAL_CONDITIONS: Final[Tuple[str, ...]] = tuple(
    [
        "early_bsl_nrem",
        "early_bsl_rem",
        "last_bsl_nrem",
        "last_bsl_rem",
        "bsl_wake",
        "bsl_rem",
        "ext_wake",
        "early_ext_wake",
        "late_ext_wake",
        "sd_wake",
        "early_sd_wake",
        "late_sd_wake",
        "nod_wake",
        "early_nod_wake",
        "late_nod_wake",
        "cow_wake",
        "early_cow_wake",
        "late_cow_wake",
        "ctn_wake",
        "early_ctn_wake",
        "late_ctn_wake",
        "early_rec_nrem",
        "early_rec_rem",
        "late_rec_nrem",
        "late_rec_rem",
        "last_rec_nrem",
        "last_rec_rem",
        "early_rec_nrem_match",
        "early_rec_rem_match",
    ]
)

NOD_STATISTICAL_CONDITIONS: Final[Tuple[str, ...]] = tuple(
    [x for x in ALL_STATISTICAL_CONDITIONS if ("cow" not in x) and ("ctn" not in x)]
)
COW_STATISTICAL_CONDITIONS: Final[Tuple[str, ...]] = tuple(
    [x for x in ALL_STATISTICAL_CONDITIONS if ("nod" not in x) and ("ctn" not in x)]
)
CTN_STATISTICAL_CONDITIONS: Final[Tuple[str, ...]] = ALL_STATISTICAL_CONDITIONS
