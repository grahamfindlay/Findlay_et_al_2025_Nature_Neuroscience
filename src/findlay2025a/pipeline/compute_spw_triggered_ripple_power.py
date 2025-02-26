import mne.filter
import numpy as np
import scipy.fftpack
import scipy.signal
import xarray as xr
from ecephys import wne

from findlay2025a import core
from findlay2025a.constants import Bands, Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")

    # Load data. ~45s
    print("Loading LFP snippets...")
    spw_lfp_snippet_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.SPW_LFP_SNIPPETS
    )
    spw_lfps = xr.load_dataarray(spw_lfp_snippet_file)

    # Get motion corrected snippets. This takes ~1m
    print("Motion correcting...")
    shifts = core.read_shifts(sglx_subject.name, experiment)
    chs = core.get_trial_shifted_channel_masks(
        shifts,
        spw_lfps,
        spw_lfps.sel(channel=spw_lfps.acronym.isin(["CA1-sp"])).channel.values,
    )
    spw_lfps = xr.concat(
        [
            spw_lfps.sel(event=evt_chans.event, channel=evt_chans).drop_vars(
                spw_lfps.channel.coords
            )
            for evt_chans in chs
        ],
        dim="event",
    )
    # To preserve channel info, at the cost of memory and speed of downstream operations: spw_lfps.where(chs, drop=True)

    # Apply ripple filter. ~4m
    # All hippocampal channels: 3h w/o parellelization. Might be faster on 2D data...
    print("Applying ripple filter...")
    spw_lfps = spw_lfps.transpose("event", "channel", "time")
    original_dtype = spw_lfps.dtype
    spw_lfps.values = mne.filter.filter_data(
        spw_lfps.values.astype(np.float64),
        spw_lfps.fs,
        Bands.RIPPLE.low,
        Bands.RIPPLE.high,
    ).astype(original_dtype)

    # Get analytic signal. ~15s
    print("Getting analytic signal...")
    n_samples = spw_lfps.time.size
    nfast = scipy.fftpack.next_fast_len(n_samples)
    spw_lfps.values = scipy.signal.hilbert(spw_lfps, N=nfast)[:, :, :n_samples]

    # Get instantaneous power.
    print("Getting instantaneous frequency and power...")
    inst_freq = xr.full_like(spw_lfps, np.nan, dtype=np.float64)
    inst_freq.values = np.unwrap(np.angle(spw_lfps))
    inst_freq = inst_freq.diff(dim="time") / (2.0 * np.pi) * inst_freq.fs
    inst_pow = xr.full_like(spw_lfps, np.nan, dtype=np.float64)
    inst_pow.values = np.log10(np.square(np.abs(spw_lfps)))

    print("Saving...")
    freq_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.SPW_TRIGGERED_RIPPLE_FREQ
    )
    inst_freq.to_netcdf(freq_file)
    power_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.SPW_TRIGGERED_RIPPLE_POWER
    )
    inst_pow.to_netcdf(power_file)

    if False:  # Useful for evaluating drift correction. So far looks flawless.
        uEv = inst_pow.mean(dim="event")
        uEv.plot.imshow(x="time", y="y")


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
