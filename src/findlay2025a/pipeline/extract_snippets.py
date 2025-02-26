import numpy as np
import xarray as xr
from ecephys import wne, xrsig

from findlay2025a import core, ripples, sharp_waves
from findlay2025a.constants import Files


def _detrend(da: xr.DataArray, blocksize=10000) -> xr.DataArray:
    """A horrible, ugly hack to detrend in blocks of 10000 trials, becuase polyfit segfaults for > ~45k trials, for reasons not obviously related to memory consumption.
    Something about numerical stability, I guess, becuase it relies on the old numpy polyfit interface."""
    blks = []
    for i in np.arange(0, da.event.size, blocksize):
        blks.append(xrsig.detrend_trialed(da.isel(event=slice(i, i + blocksize))))
    return xr.concat(blks, dim="event")


def do_experiment(
    sglx_subject: wne.sglx.SGLXSubject,
    experiment: str,
):
    nb = core.get_project("seahorse")

    lfp = core.open_hippocampal_lfps(sglx_subject.name, experiment)
    spws = sharp_waves.read_spws(sglx_subject.name, experiment, kind="postprocessed")
    rips = ripples.read_ripples(sglx_subject.name, experiment, kind="postprocessed")

    print("Loading LFP")
    lfp = lfp.load()
    print("Loading KCSD")
    kcsd = core.open_kcsd(sglx_subject.name, experiment).load()
    kcsd.attrs["fs"] = lfp.fs

    print("Extracting SPW LFP snippets")
    spw_lfps, in_spw_lfps = xrsig.make_trialed(
        lfp, pre=0.350, post=0.350, event_times=spws.pk_time.values
    )
    if not in_spw_lfps.all():
        print("Warning: Not all events are in bounds")

    spw_lfps = _detrend(spw_lfps)
    spw_lfp_snippet_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.SPW_LFP_SNIPPETS
    )
    spw_lfps.to_netcdf(spw_lfp_snippet_file)

    print("Extracting ripple LFP snippets")
    rip_lfps, in_rip_lfps = xrsig.make_trialed(
        lfp, pre=0.350, post=0.350, event_times=rips.pk_time.values
    )
    if not in_rip_lfps.all():
        print("Warning: Not all events are in bounds")

    rip_lfps = _detrend(rip_lfps)
    rip_lfp_snippet_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_LFP_SNIPPETS
    )
    rip_lfps.to_netcdf(rip_lfp_snippet_file)

    print("Extracting ripple KCSD snippets")
    rip_csds, in_rip_csds = xrsig.make_trialed(
        kcsd, pre=0.350, post=0.350, event_times=rips.pk_time.values
    )
    if not in_rip_csds.all():
        print("Warning: Not all events are in bounds")

    rip_csd_snippet_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_CSD_SNIPPETS
    )
    rip_csds.to_netcdf(rip_csd_snippet_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject.name} {experiment}")
        do_experiment(sglx_subject, experiment)
