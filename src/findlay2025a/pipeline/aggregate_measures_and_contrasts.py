import pandas as pd

from findlay2025a import core
from findlay2025a.constants import Files


def do_project():
    nb = core.get_project("seahorse")

    cpc = pd.read_parquet(nb.get_project_file(Files.CX_BANDPOWER_CONTRASTS))
    cpm = pd.read_parquet(nb.get_project_file(Files.CX_BANDPOWER_MEANS))

    hpm = pd.read_parquet(nb.get_project_file(Files.HIPPOCAMPAL_BANDPOWER_MEANS))
    hps = pd.read_parquet(nb.get_project_file(Files.HIPPOCAMPAL_BANDPOWER_SUM))

    csc = pd.read_parquet(nb.get_project_file(Files.SPW_CONDITION_CONTRASTS))
    csr = pd.read_parquet(nb.get_project_file(Files.SPW_CONDITION_RATES))
    csm = pd.read_parquet(nb.get_project_file(Files.SPW_CONDITION_MEANS))

    dsc = pd.read_parquet(nb.get_project_file(Files.DSPK_CONDITION_CONTRASTS))
    dsr = pd.read_parquet(nb.get_project_file(Files.DSPK_CONDITION_RATES)).drop(
        columns=["duration"]
    )  # Already present in SPW condition rates
    dsm = pd.read_parquet(nb.get_project_file(Files.DSPK_CONDITION_MEANS))

    rpc = pd.read_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_CONTRASTS))
    rpr = pd.read_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_RATES)).drop(
        columns=["duration"]
    )  # Already present in SPW condition rates
    rpm = pd.read_parquet(nb.get_project_file(Files.RIPPLE_CONDITION_MEANS))

    emm = pd.read_parquet(nb.get_project_file(Files.EMG_CONDITION_MEANS))
    ems = pd.read_parquet(nb.get_project_file(Files.EMG_CONDITION_SUMS))

    cxspm = pd.read_parquet(nb.get_project_file(Files.CORTICAL_SPINDLE_CONDITION_MEANS))
    cxspr = pd.read_parquet(
        nb.get_project_file(Files.CORTICAL_SPINDLE_CONDITION_RATES)
    ).drop(columns=["duration"])  # Already present in SPW condition rates

    measures = pd.concat(
        [cpm, hpm, hps, csm, csr, dsm, dsr, rpm, rpr, emm, ems, cxspr, cxspm],
        axis="columns",
    )
    measures.to_parquet(nb.get_project_file(Files.COMBINED_CONDITION_MEASURES))

    contrasts = pd.concat([cpc, csc, dsc, rpc], axis="columns")
    contrasts.to_parquet(nb.get_project_file(Files.COMBINED_CONDITION_CONTRASTS))
