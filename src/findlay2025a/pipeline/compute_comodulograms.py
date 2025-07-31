import logging

import numpy as np
import pactools
import xarray as xr
from tqdm.auto import tqdm

from ecephys import wne
from findlay2025a import core, hypnograms
from findlay2025a.constants import Files

logger = logging.getLogger(__name__)

DURATION = 600
N_SURROGATES = 200
LOW_FQ_RANGE = np.linspace(1, 10, 50)
LOW_FQ_WIDTH = 3.2


def get_pac_lfp(
    experiment: str, sglx_subject: wne.sglx.SGLXSubject, state: str, duration: float
) -> xr.DataArray:
    full_hg = hypnograms.load_consolidated_hypnogram(
        experiment, sglx_subject.name, simplify=True, clean=True
    )
    hg = full_hg.keep_states([state]).keep_first(duration)
    lf = core.open_hippocampal_lfps(sglx_subject.name, experiment)
    return xr.concat(
        [
            lf.sel(time=slice(bout.start_time, bout.end_time))
            for bout in hg.itertuples()
        ],
        dim="time",
    )


def comodulogram_to_dataarray(c: pactools.Comodulogram) -> xr.DataArray:
    return xr.DataArray(
        c.comod_,
        dims=("driver_frequency", "signal_frequency"),
        coords={
            "driver_frequency": c.low_fq_range,
            "signal_frequency": c.high_fq_range,
        },
    )


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str, state: str):
    assert N_SURROGATES > 1, "Must have at least 2 surrogates, preferably 200"
    nb = core.get_project("seahorse")
    lf = get_pac_lfp(experiment, sglx_subject, state, DURATION)

    cdict = {}
    for chan in tqdm(lf.channel.values, desc="Channel..."):
        lf_chan = lf.sel(channel=chan)
        cdict[chan] = pactools.Comodulogram(
            fs=lf_chan.fs,
            low_fq_range=LOW_FQ_RANGE,
            low_fq_width=LOW_FQ_WIDTH,
            method="duprelatour",
            progress_bar=False,
            n_surrogates=N_SURROGATES,
            n_jobs=-1,
        )
        cdict[chan].fit(lf_chan.values)

    comods = xr.concat(
        [
            comodulogram_to_dataarray(cdict[chan])
            .assign_coords(**lf.sel(channel=chan)["channel"].coords)
            .expand_dims("channel")
            for chan in cdict
        ],
        dim="channel",
    )
    comods.to_netcdf(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, Files.COMODULOMGRAMS(state)
        )
    )
    surrogate_maxes = xr.concat(
        [
            xr.DataArray(cdict[chan].surrogate_max_, dims="surrogate")
            .assign_coords(channel=chan)
            .expand_dims("channel")
            for chan in cdict
        ],
        dim="channel",
    )
    surrogate_maxes.to_netcdf(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, Files.COMODULOMGRAM_SURROGATE_MAXES(state)
        )
    )


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment, "NREM")
        do_experiment(sglx_subject, experiment, "REM")
