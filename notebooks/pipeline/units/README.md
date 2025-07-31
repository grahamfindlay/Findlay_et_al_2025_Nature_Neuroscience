Run notebooks in the following order:
1. `compute_acgs.ipynb`
2. `compute_acg_metrics.ipynb`
3. `agg_cell_metrics_assign_unit_quality.ipynb`: 
Produces `mps_metrics.pqt` (`mps` = "multiprobe sorting"), `aggregated_cell_metrics.pqt` and `cluster_quality.pqt`
4. `assign_cell_types.ipynb`
Produces `cell_types.pqt`, updates `aggregated_cell_metrics.pqt`: 
5. `compute_peths.ipynb`. 
6. `compute_frs_by_condition.ipynb`: 
Produces `firing_rates_by_condition.pqt`, which includes clusters of all quality (including multiunit activity). 
Also produces `frs.with_meta.sua_moderate.pqt`, which only includes single unit clusters. 
7. `estimate_on_frs_by_condition.ipynb`: 
Removes detected OFF periods to obtain adjusted ON firing rates, produces `estm_on_frs.with_meta.sua_moderate.pqt`.