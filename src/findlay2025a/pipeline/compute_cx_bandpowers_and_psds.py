# TODO: Use just a subset of channels, e.g. some ROI chans.
from typing import Final

import xarray as xr
from ecephys import npsig, wne, xrsig

from findlay2025a import core, hypnograms
from findlay2025a.constants import Bands

STFT_CHUNK_DURATION: Final[float] = 4  # Seconds
DO_BANDS: Final[list[Bands]] = [
    Bands.DELTA,
    Bands.THETA,
    Bands.SIGMA,
    Bands.GAMMA,
]
Q: Final[int] = 2  # Temporal decimation factor


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    # Reference cortical LFPs to white matter
    cx_lfps = core.open_cortical_lfps(
        sglx_subject.name, experiment, drop_bad_channels=True
    )
    wm_lfps = core.open_white_matter_lfps(
        sglx_subject.name, experiment, drop_bad_channels=True
    )
    lfps = cx_lfps - wm_lfps.mean(dim="channel")
    lfps.attrs = cx_lfps.attrs

    lfps = xrsig.antialiasing_filter(lfps, q=Q)
    lfps = lfps.isel(time=slice(None, None, Q))
    lfps.attrs["fs"] /= Q

    _, (n_fft, _) = npsig.get_n_fft(lfps.fs, s=STFT_CHUNK_DURATION)
    spgs = xrsig.stft_psd(lfps, n_fft=n_fft)

    for band in DO_BANDS:
        pwr = (
            spgs.sel(frequency=slice(*band.value))
            .sum(dim="frequency")
            .mean(dim="channel")
        )
        savefile = core.get_cortical_bandpower_file(sglx_subject.name, experiment, band)
        pwr.to_netcdf(savefile)

    hgs = hypnograms.load_statistical_condition_hypnograms(
        experiment, sglx_subject.name
    )
    psds = xr.Dataset(
        {c: spgs.sel(time=hgs[c].covers_time(spgs.time)).mean(dim="time") for c in hgs}
    )
    savefile = core.get_cortical_psds_file(sglx_subject.name, experiment)
    psds.to_netcdf(savefile)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
