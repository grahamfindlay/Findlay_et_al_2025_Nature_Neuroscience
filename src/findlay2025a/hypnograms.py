import warnings

from ecephys import hypnogram, wne
from findlay2025a import constants, core
from wisc_ecephys_tools.rats import cnd_hgs


def load_statistical_condition_hypnograms(
    experiment: str,
    subject: str,
    include_full_conservative: bool = False,
    use_all_new: bool = False,
) -> dict[str, hypnogram.FloatHypnogram]:
    nb = core.get_project("seahorse")
    probe = core.get_spw_probe(experiment, subject)
    new_hgs = cnd_hgs.load_statistical_condition_hypnograms(
        nb.get_experiment_subject_file(
            experiment,
            subject,
            f"{probe}.condition_hypnograms.parquet",
        )
    )
    conditions = (
        constants.CORE_STATISTICAL_CONDITIONS + ["full_conservative"]
        if include_full_conservative
        else constants.CORE_STATISTICAL_CONDITIONS
    )  # TODO: Replace `CORE_...` with expected conditions based on `experiment`.
    hgs = {}
    for c in conditions:
        if use_all_new or (c in ["early_ext", "late_ext"]):
            hgs[c] = new_hgs[c]
        else:
            try:
                hgs[c] = hypnogram.FloatHypnogram.from_htsv(
                    nb.get_experiment_subject_file(
                        experiment, subject, f"{c}_hypnogram.htsv"
                    )
                )
            except FileNotFoundError:
                pass
    return hgs


def load_consolidated_hypnogram(
    experiment: str, subject: str, simplify: bool = True, clean: bool = True
) -> hypnogram.FloatHypnogram:
    s3 = core.get_project("shared")
    opts = s3.load_experiment_subject_params(experiment, subject)
    try:
        hg = wne.utils.load_consolidated_hypnogram(
            s3,
            experiment,
            subject,
            opts["hypnogram_probe"],
            simplify=simplify,
            clean=clean,
        )
    except FileNotFoundError:
        warnings.warn(
            f"Consolidated hypnogram for {experiment} {subject} {opts['hypnogram_probe']} "
            "not found. Please generate this file. The use of probe-agnostic hypnograms "
            "is deprecated. Trying to load one anyways.",
            DeprecationWarning,
            stacklevel=2,
        )
        f = s3.get_experiment_subject_file(experiment, subject, wne.Files.HYPNOGRAM)
        hg = hypnogram.FloatHypnogram.from_htsv(f)
        if simplify:
            hg = hg.replace_states(wne.SIMPLIFIED_STATES)
            if clean:
                # Cleaning is already done, but becuase we are simplifying states, it might need
                # to be done again. It will change NaNs to NoData. Is that expected downstream somewhere?
                hg = hypnogram.FloatHypnogram.clean(hg._df)

    if clean:
        warnings.warn(
            "The use of `clean` is deprecated. Please confirm that it is not needed, "
            "then set it to `False`.",
            DeprecationWarning,
            stacklevel=2,
        )

    return hg
