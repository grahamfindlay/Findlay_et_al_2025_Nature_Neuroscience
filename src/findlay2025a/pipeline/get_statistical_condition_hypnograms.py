import pandas as pd

import wisc_ecephys_tools as wet
from ecephys import wne
from findlay2025a import core
from wisc_ecephys_tools.rats import cnd_hgs, exp_hgs


def do_probe(subject: wne.sglx.SGLXSubject, experiment: str, probe: str):
    s3 = core.get_project("shared")
    assert probe == core.get_spw_probe(experiment, subject.name)

    # IMPT: We will use scoring from periods where we may not have actually been able to collect SPWs.
    # We are very careful to handle this later, e.g. when getting event rates!
    lbrl_hg = exp_hgs.get_liberal_hypnogram(s3, experiment, subject, probe)

    # Load a more conservative hypnogram with artifacts and such that we might want to know about.
    cons_hg = exp_hgs.load_hypnogram(
        s3,
        experiment,
        subject,
        probe,
        include_ephyviewer_edits=True,
        include_sorting_nodata=False,
        include_lf_consolidated_artifacts=True,
        include_ap_consolidated_artifacts=False,
        include_lf_sglx_filetable_nodata=True,
        include_ap_sglx_filetable_nodata=False,
        simplify=True,
        fallback=True,
    )  # TODO: Would be nice to use all sources. Must check impact of this on results.

    return cnd_hgs.compute_statistical_condition_hypnograms(
        lbrl_hg, cons_hg, experiment, subject
    )


def do_experiment_subject(
    sglx_subject: wne.sglx.SGLXSubject,
    experiment: str,
    probes: tuple[str, ...],
    verbose: bool = True,
    save: bool = True,
):
    prb_hgs = {prb: do_probe(sglx_subject, experiment, prb) for prb in probes}

    if save:
        nb = wet.get_sglx_project("seahorse")
        for prb, hgs in prb_hgs.items():
            fpath = nb.get_experiment_subject_file(
                experiment,
                sglx_subject.name,
                f"{prb}.condition_hypnograms.parquet",
            )
            cnd_hgs.save_statistical_condition_hypnograms(hgs, fpath)
    if len(prb_hgs) < 2:
        return None, None, prb_hgs

    consensus_hgs, consensus_df = cnd_hgs.get_consensus(prb_hgs)
    if verbose:
        pd.set_option("display.max_rows", 100)
        print(consensus_df)
        pd.reset_option("display.max_rows")
    if save:
        fpath = nb.get_experiment_subject_file(
            experiment,
            sglx_subject.name,
            "consensus_condition_hypnograms.parquet",
        )
        cnd_hgs.save_statistical_condition_hypnograms(consensus_hgs, fpath)
    return consensus_hgs, consensus_df, prb_hgs


def do_project(verbose: bool = True, save: bool = True):
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        probes = (core.get_spw_probe(experiment, sglx_subject.name),)
        print(f"Processing {sglx_subject.name} {experiment}, {probes}")
        do_experiment_subject(sglx_subject, experiment, probes, verbose, save)


if __name__ == "__main__":
    do_project()
