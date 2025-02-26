import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from ecephys import wne, xrsig
from mpl_toolkits.axes_grid1 import ImageGrid

from findlay2025a import core
from findlay2025a.constants import Files


def plot_laminar_comodulograms(
    comods: xr.DataArray,
    surrogate_maxes: xr.DataArray = None,
    p_value=0.05,
    ncols=10,
    size=4,
    aspect=0.8,
    ax_title_kwargs=dict(fontsize=10, color="w", y=0.9, horizontalalignment="right"),
    share_vlim=False,
):
    if share_vlim:
        (vmin, vmax) = (float(comods.min()), float(comods.max()))
    else:
        (vmin, vmax) = (None, None)

    nrows = int(np.ceil(len(comods.channel) / ncols))
    fig = plt.figure(figsize=(size * aspect * ncols, size * nrows))
    grid = ImageGrid(fig, 111, (nrows, ncols), aspect=False)
    for chan, y, ax in zip(comods.channel.values[::-1], comods.y.values[::-1], grid):
        c = comods.sel(channel=chan)
        extent = (
            c.driver_frequency.min(),
            c.driver_frequency.max(),
            c.signal_frequency.min(),
            c.signal_frequency.max(),
        )
        ax.imshow(
            c.values.T,
            origin="lower",
            extent=extent,
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
            interpolation="none",
        )
        if surrogate_maxes is not None:
            sm = surrogate_maxes.sel(channel=chan)
            percentile = 100 * (1 - p_value)
            level = np.percentile(sm, percentile)
            ax.contour(
                c.driver_frequency.values,
                c.signal_frequency.values,
                c.values.T,
                levels=[level],
                colors="w",
                origin="lower",
            )
        ax.set_title(f"LF{chan}: {y}um", **ax_title_kwargs)
    return fig


def do_experiment(
    sglx_subject: wne.sglx.SGLXSubject,
    experiment: str,
    state: str,
    do_extras: bool = False,
):
    nb = core.get_project("seahorse")
    comods = xr.load_dataarray(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, Files.COMODULOMGRAMS(state)
        )
    )
    surrogate_maxes = xr.load_dataarray(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, Files.COMODULOMGRAM_SURROGATE_MAXES(state)
        )
    )

    fig = plot_laminar_comodulograms(comods, surrogate_maxes=surrogate_maxes)
    fig.savefig(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, f"plots/{state}_comodulograms.png"
        )
    )
    plt.close(fig)

    fig = plot_laminar_comodulograms(comods, share_vlim=True)
    fig.savefig(
        nb.get_experiment_subject_file(
            experiment,
            sglx_subject.name,
            f"plots/{state}_comodulograms_shared_vlim.png",
        )
    )
    plt.close(fig)

    if not do_extras:
        return

    modal_driver = float(
        comods.driver_frequency[
            comods.max(dim="signal_frequency").argmax(dim="driver_frequency")
        ].median()
    )
    mods = comods.sel(driver_frequency=modal_driver, method="nearest")
    mods_norm = mods / mods.max(dim="signal_frequency")

    xrsig.plot_laminar_image_horizontal(mods_norm)
    fig = plt.gcf()
    fig.savefig(
        nb.get_experiment_subject_file(
            experiment, sglx_subject.name, f"plots/{state}_laminar_driver_pac.png"
        ),
        bbox_inches="tight",
    )
    plt.close(fig)

    xrsig.plot_laminar_image_horizontal(mods)
    fig = plt.gcf()
    fig.savefig(
        nb.get_experiment_subject_file(
            experiment,
            sglx_subject.name,
            f"plots/{state}_laminar_driver_pac_shared_vlim.png",
        ),
        bbox_inches="tight",
    )
    plt.close(fig)


def do_project():
    for sglx_subject, experiment in core.yield_sglx_subject_experiment_pairs():
        print(f"Doing {sglx_subject} {experiment}")
        do_experiment(sglx_subject, experiment, "NREM")
        do_experiment(sglx_subject, experiment, "REM")
