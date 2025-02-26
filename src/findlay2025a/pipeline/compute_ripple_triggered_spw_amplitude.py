import xarray as xr
from ecephys import wne

from findlay2025a import core
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")

    # Load data. ~30s
    print("Loading CSD snippets...")
    rip_csd_snippet_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_CSD_SNIPPETS
    )
    rip_csds = xr.load_dataarray(rip_csd_snippet_file)

    # We need to take acronym info from LFPs
    rip_lfp_snippet_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_LFP_SNIPPETS
    )
    rip_lfps = xr.open_dataarray(rip_lfp_snippet_file)

    rip_lfps.channel.load()
    common_coords = rip_lfps.sel(channel=rip_csds.channel).channel
    assert (common_coords.channel.size == rip_csds.channel.size) and (
        common_coords.channel.size <= rip_lfps.channel.size
    )
    rip_csds = rip_csds.assign_coords(channel=common_coords)

    # Get motion corrected snippets. ~1m
    print("Motion correcting...")
    shifts = core.read_shifts(sglx_subject.name, experiment)
    chs = core.get_trial_shifted_channel_masks(
        shifts,
        rip_csds,
        rip_csds.sel(channel=rip_csds.acronym.isin(["CA1-sr"])).channel.values,
    )
    rip_csds = rip_csds.where(chs, drop=True)

    print("Saving...")
    sink_amps = -rip_csds.min(dim="channel")
    amp_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_TRIGGERED_SPW_AMP
    )
    sink_amps.to_netcdf(amp_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
