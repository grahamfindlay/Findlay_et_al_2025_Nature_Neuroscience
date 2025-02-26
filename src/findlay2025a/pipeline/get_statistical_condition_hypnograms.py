from ecephys import wne

from findlay2025a import core, hypnograms


def get_statistical_condition_hypnograms(
    sglx_subject: wne.sglx.SGLXSubject, experiment: str
):
    s3 = core.get_project("shared")

    # IMPT: We will use scoring from periods where we may not have actually been able to collect SPWs.
    # We are very careful to handle this later, e.g. when getting event rates!
    exp_hg = wne.sglx.utils.load_reconciled_float_hypnogram(
        s3,
        experiment,
        sglx_subject,
        probes=[],
        sources=[],
        reconcile_ephyviewer_edits=True,
        simplify=True,
    )

    # Load a more conservative hypnogram with artifacts and such that we might want to know about.
    prb_hg = wne.sglx.utils.load_reconciled_float_hypnogram(
        s3,
        experiment,
        sglx_subject,
        probes=[core.get_spw_probe(experiment, sglx_subject.name)],
        sources=[
            "lf"
        ],  # TODO: Would be nice to use "sorting" source too, but currently not supported.
        reconcile_ephyviewer_edits=True,
        simplify=True,
        alias="full",
        sorting="sorting",
    )

    sc_hgs = hypnograms.compute_statistical_condition_hypnograms(
        exp_hg, prb_hg, experiment, sglx_subject
    )

    return sc_hgs


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    sc_hgs = get_statistical_condition_hypnograms(sglx_subject, experiment)

    for c, c_hg in sc_hgs.items():
        c_hg.write_htsv(
            nb.get_experiment_subject_file(
                experiment, sglx_subject.name, f"{c}_hypnogram.htsv"
            )
        )


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Processing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
