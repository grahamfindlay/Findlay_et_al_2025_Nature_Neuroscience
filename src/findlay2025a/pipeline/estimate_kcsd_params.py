import numpy as np
from ecephys import wne, xrsig

from findlay2025a import core
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    """Requires the following fields to be defined in experiment_params.json...
    spwrProbe: Use this probe for parameter estimation
    probes->imecX->badChannels: List of bad channels.
    probes->imexX->structureBounds->hippocampus: Borders of hippocampus, in um
    """
    # Get LFP
    t1, t2 = core.get_estimation_bounds(sglx_subject, experiment)
    lf = core.open_hippocampal_lfps(sglx_subject.name, experiment).sel(
        time=slice(t1, t2)
    )

    # Estimate CSD
    s3 = core.get_project("shared")
    params = s3.load_experiment_subject_params(experiment, sglx_subject.name)
    prb = params["spwrProbe"]
    bad_chans = np.array(params["probes"][prb]["badChannels"])
    csd = xrsig.kernel_current_source_density(lf, drop=bad_chans, do_lcurve=True)

    # Store CSD params
    nb = core.get_project("seahorse")
    csd_params_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.KCSD_PARAMS
    )
    np.savez(
        csd_params_file,
        estm_start_time=t1,
        estm_end_time=t2,
        electrode_pitch=csd.pitch_mm,
        xmin=csd.kcsd.xmin,
        xmax=csd.kcsd.xmax,
        n_estm=csd.kcsd.n_estm,
        gdx=csd.kcsd.gdx,
        lambd=csd.kcsd.lambd,
        R=csd.kcsd.R,
        channels_used=csd.channel.values,
        channels_omitted=bad_chans,
        ele_pos=csd.kcsd.ele_pos,
    )


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
