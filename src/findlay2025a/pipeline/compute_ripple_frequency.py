import mne.filter
import numpy as np
import scipy.fftpack
import scipy.signal
import xarray as xr
from ecephys import wne

from findlay2025a import core, ripples
from findlay2025a.constants import Bands, Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")

    # Load data. ~2m
    print("Loading LFP snippets...")
    rip_lfp_snippet_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_LFP_SNIPPETS
    )
    rip_lfps = xr.load_dataarray(rip_lfp_snippet_file)

    # Get motion corrected snippets. ~1.5m
    print("Motion correcting...")
    rips = ripples.read_ripples(sglx_subject.name, experiment, kind="postprocessed")
    ripple_pk_chans = (
        rips.set_index("pk_time")["pk_chan_id"].loc[rip_lfps.event.values].reset_index()
    )
    rip_lfps = xr.concat(
        [
            rip_lfps.sel(event=e.pk_time, channel=int(e.pk_chan_id))
            for e in ripple_pk_chans.itertuples()
        ],
        dim="event",
    )
    # Discards channel dimension, but keeps channel as a coord on event dimm

    # Apply ripple filter. ~11s
    print("Applying ripple filter...")
    original_dtype = rip_lfps.dtype
    rip_lfps.values = mne.filter.filter_data(
        rip_lfps.values.astype(np.float64),
        rip_lfps.fs,
        Bands.RIPPLE.low,
        Bands.RIPPLE.high,
    ).astype(original_dtype)

    # Get analytic signal. ~1s
    print("Getting analytic signal...")
    n_samples = rip_lfps.time.size
    nfast = scipy.fftpack.next_fast_len(n_samples)
    rip_lfps.values = scipy.signal.hilbert(rip_lfps, N=nfast)[
        :, :n_samples
    ]  # get analytic signal. 4m15s

    print("Getting instantaneous frequency and power...")
    inst_freq = xr.full_like(rip_lfps, np.nan, dtype=np.float64)
    inst_freq.values = np.unwrap(np.angle(rip_lfps))
    inst_freq = inst_freq.diff(dim="time") / (2.0 * np.pi) * inst_freq.fs
    inst_pow = xr.full_like(rip_lfps, np.nan, dtype=np.float64)
    inst_pow.values = np.log10(np.square(np.abs(rip_lfps)))

    print("Saving...")
    freq_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_FREQ
    )
    inst_freq.to_netcdf(freq_file)
    power_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.RIPPLE_POWER
    )
    inst_pow.to_netcdf(power_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
