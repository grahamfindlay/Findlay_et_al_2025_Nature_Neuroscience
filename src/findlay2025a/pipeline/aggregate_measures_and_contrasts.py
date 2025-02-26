import pandas as pd

from findlay2025a import core
from findlay2025a.constants import Files


def do_project():
    nb = core.get_project("seahorse")

    cpc = pd.read_parquet(nb.get_project_file(Files.CX_BANDPOWER_CONTRASTS))
    cpm = pd.read_parquet(nb.get_project_file(Files.CX_BANDPOWER_MEANS))
    cps = pd.read_parquet(nb.get_project_file(Files.CX_BANDPOWER_SUM))

    hpc = pd.read_parquet(nb.get_project_file(Files.HIPPOCAMPAL_BANDPOWER_CONTRASTS))
    hpm = pd.read_parquet(nb.get_project_file(Files.HIPPOCAMPAL_BANDPOWER_MEANS))
    hps = pd.read_parquet(nb.get_project_file(Files.HIPPOCAMPAL_BANDPOWER_SUM))

    csc = pd.read_parquet(nb.get_project_file(Files.SPW_CONDITION_CONTRASTS))
    csr = pd.read_parquet(nb.get_project_file(Files.SPW_CONDITION_RATES))
    csm = pd.read_parquet(nb.get_project_file(Files.SPW_CONDITION_MEANS))
    css = pd.read_parquet(nb.get_project_file(Files.SPW_CONDITION_SUMS))

    dsc = pd.read_parquet(nb.get_project_file(Files.DSPK_CONDITION_CONTRASTS))
    dsr = pd.read_parquet(nb.get_project_file(Files.DSPK_CONDITION_RATES)).drop(
        columns=["duration"]
    )  # Already present in SPW condition rates
    dsm = pd.read_parquet(nb.get_project_file(Files.DSPK_CONDITION_MEANS))
    dss = pd.read_parquet(nb.get_project_file(Files.DSPK_CONDITION_SUMS))

    rpc = pd.read_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_CONTRASTS))
    rpr = pd.read_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_RATES)).drop(
        columns=["duration"]
    )  # Already present in SPW condition rates
    rpm = pd.read_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_MEANS))
    rps = pd.read_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_SUMS))

    emc = pd.read_parquet(nb.get_project_file(Files.EMG_CONDITION_CONTRASTS))
    emm = pd.read_parquet(nb.get_project_file(Files.EMG_CONDITION_MEANS))
    ems = pd.read_parquet(nb.get_project_file(Files.EMG_CONDITION_SUMS))

    measures = pd.concat(
        [cpm, cps, hpm, hps, csm, css, csr, dsm, dss, dsr, rpm, rps, rpr, emm, ems],
        axis="columns",
    )
    measures.to_parquet(nb.get_project_file(Files.COMBINED_CONDITION_MEASURES))

    contrasts = pd.concat([cpc, hpc, csc, dsc, rpc, emc], axis="columns")
    contrasts.to_parquet(nb.get_project_file(Files.COMBINED_CONDITION_CONTRASTS))
