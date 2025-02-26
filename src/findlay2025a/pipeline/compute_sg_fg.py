import numpy as np
from ecephys import utils, wne, xrsig

from findlay2025a import core, hypnograms
from findlay2025a.constants import Bands


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")
    s3 = core.get_project("shared")
    lf = core.open_hippocampal_lfps(sglx_subject.name, experiment)

    # Get ROIs
    y = lf.y.values
    channel_info = lf.channel.to_dataframe()
    median_depths = channel_info.groupby("acronym")["y"].median()
    roi_labels = ["CA1-so", "CA1-sp", "CA1-sr", "CA1-slm", "DG-dl-ml"]
    roi_depths = [
        y[utils.find_nearest(y, median_depths.loc[roi])] for roi in roi_labels
    ]
    opts = s3.load_experiment_subject_params(experiment, sglx_subject.name)
    if "ventral_dentate_source_hi" in opts:
        ca3_depth = (
            opts["ventral_dentate_source_hi"] + opts["dorsal_dentate_source_lo"]
        ) / 2
    else:
        ca3_depth = (
            opts["dorsal_dentate_source_hi"] + opts["dorsal_dentate_source_lo"]
        ) / 2
    ca3_depth = y[utils.find_nearest(y, ca3_depth)]
    roi_depths.append(ca3_depth)
    roi_labels.append("CA3")
    roi_chans = channel_info.set_index("y").loc[roi_depths].channel.to_list()

    # Drop artifacts
    hgs = hypnograms.load_statistical_condition_hypnograms(
        experiment, sglx_subject.name, include_full_conservative=True
    )
    hg = hgs.pop("full_conservative").drop_states(core.ARTIFACT_STATES + ["NoData"])
    lf = lf.isel(time=hg.covers_time(lf["time"]))

    # Load data into memory
    print("Loading data...")
    lf_rois = (
        lf.sel(channel=roi_chans).assign_coords({"roi": ("channel", roi_labels)}).load()
    )

    # Compute STFT
    print("Computing STFT...")
    stft = xrsig.complex_stft(lf_rois, n_fft=int(lf.fs), hop_len=int(lf.fs * 0.5))

    # Get power
    print("Computing power...")
    spgs = np.abs(stft) ** 2

    # Get slow and fast gamma, with ratio
    print("Saving slow and fast gamma...")
    sg = spgs.sel(frequency=slice(*Bands.AERY_JONES_SLOW_GAMMA)).sum(dim="frequency")
    fg = spgs.sel(frequency=slice(*Bands.AERY_JONES_FAST_GAMMA)).sum(dim="frequency")
    sg.to_netcdf(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, "AeryJones_slow_gamma.nc"
        )
    )
    fg.to_netcdf(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, "AeryJones_fast_gamma.nc"
        )
    )


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
