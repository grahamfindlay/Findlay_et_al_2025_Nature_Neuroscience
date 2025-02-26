import pandas as pd
import scipy.stats
import xarray as xr
from ecephys import wne

from findlay2025a import core, ripples, sharp_waves
from findlay2025a.constants import Files


def do_experiment(sglx_subject: wne.sglx.SGLXSubject, experiment: str):
    nb = core.get_project("seahorse")

    # Get ripple-triggered properties
    rips = ripples.read_ripples(sglx_subject.name, experiment, kind="postprocessed")
    n_rips = len(rips)
    rtr_inst_freq = xr.load_dataarray(
        nb.get_experiment_subject_file(experiment, sglx_subject.name, Files.RIPPLE_FREQ)
    )
    rtr_inst_pow = xr.load_dataarray(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, Files.RIPPLE_POWER
        )
    )
    rts_sink_amps = xr.load_dataarray(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, Files.RIPPLE_TRIGGERED_SPW_AMP
        )
    )

    df1 = (
        rtr_inst_pow.sel(time=slice(-0.05, 0.05))
        .mean(dim="time")
        .to_dataframe(name="dB")[["dB"]]
    )
    df1["zdB"] = scipy.stats.zscore(df1["dB"])
    df2 = (
        rtr_inst_freq.sel(time=slice(-0.05, 0.05))
        .median(dim="time")
        .to_dataframe(name="freq")[["freq"]]
    )
    df2["zfreq"] = scipy.stats.zscore(df2["freq"])
    df3 = (
        rts_sink_amps.sel(time=slice(-0.05, 0.05))
        .mean(dim="time")
        .to_dataframe(name="sink")[["sink"]]
    )
    df3["zsink"] = scipy.stats.zscore(df3["sink"])

    df123 = pd.concat([df1, df2, df3], axis=1)
    rips = rips.drop(columns=df123.columns, errors="ignore")
    rips = rips.merge(
        df123.reset_index().rename(columns={"event": "pk_time"}),
        on="pk_time",
        how="outer",
    )
    assert len(rips) == n_rips, "Ripples were lost"
    rip_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.POSTPROCESSED_RIPPLES
    )
    rips.to_parquet(rip_file)

    # Get SPW-triggered properties
    spws = sharp_waves.read_spws(sglx_subject.name, experiment, kind="postprocessed")
    n_spws = len(spws)
    str_inst_freq = xr.load_dataarray(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, Files.SPW_TRIGGERED_RIPPLE_FREQ
        )
    )
    str_inst_pow = xr.load_dataarray(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, Files.SPW_TRIGGERED_RIPPLE_POWER
        )
    )

    df4 = (
        str_inst_freq.sel(time=slice(-0.05, 0.05))
        .median(dim="time")
        .mean(dim="channel")
        .to_dataframe(name="ripple_freq")[["ripple_freq"]]
    )
    df4["ripple_zfreq"] = scipy.stats.zscore(df4["ripple_freq"])
    df5 = (
        str_inst_pow.sel(time=slice(-0.05, 0.05))
        .median(dim="time")
        .mean(dim="channel")
        .to_dataframe(name="ripple_dB")[["ripple_dB"]]
    )
    df5["ripple_zdB"] = scipy.stats.zscore(df5["ripple_dB"])

    df45 = pd.concat([df4, df5], axis=1)
    spws = spws.drop(columns=df45.columns, errors="ignore")
    spws = spws.merge(
        df45.reset_index().rename(columns={"event": "pk_time"}),
        on="pk_time",
        how="outer",
    )
    assert len(spws) == n_spws, "SPWs were lost"
    spw_file = nb.get_experiment_subject_file(
        experiment, sglx_subject.name, Files.POSTPROCESSED_SPWS
    )
    spws.to_parquet(spw_file)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment)
