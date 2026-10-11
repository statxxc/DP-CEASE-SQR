# CEASE-SQR reproducibility code

This archive accompanies the manuscript on communication-efficient distributed
smoothed quantile regression and its differentially private extension.

## Files

- `simulation_methods.R`: standalone code for all nonprivate and private
  Monte Carlo experiments and the associated ggplot2 figures.
- `real_data_application.R`: standalone code for the 2024 ACS California PUMS
  application, including cohort construction, descriptive statistics,
  nonprivate analysis, private analysis, tables, and figures.
- `README.md`: these reproduction instructions.

Neither R file calls `source()` or depends on another code file in the archive.

## Software

The code was checked with R 4.5.3 on Windows. Install the two non-base packages:

```r
install.packages(c("ggplot2", "patchwork"))
```

The simulation also uses the base/recommended `parallel`, `stats`, and
`grDevices` facilities.

## Simulation study

A quick five-replication check using one worker process is:

```text
Rscript --vanilla simulation_methods.R both quick 1 5
```

The manuscript setting of 200 Monte Carlo replications is:

```text
Rscript --vanilla simulation_methods.R both full 6 200
```

The two components may be run separately:

```text
Rscript --vanilla simulation_methods.R nonprivate full 6 200
Rscript --vanilla simulation_methods.R private full 6 200
```

To redraw figures from completed full-summary CSV files without rerunning the
Monte Carlo experiments:

```text
Rscript --vanilla simulation_methods.R both plot_only 1
```

Simulation results are written to `simulation_results/`. The principal
manuscript figures are:

- `bandwidth_compare_full_nonprivate.pdf`
- `bandwidth_compare_full_private.pdf`

The nonprivate study compares CEASE-SQR (h,h), CEASE-SQR (h,b), CSL-SQR (h,h),
CSL-SQR (h,b), CEASE-LR, and CSL-LR. The private study compares the clean and
DP versions of CEASE-SQR under one- and two-bandwidth specifications. Full mode
uses 200 fresh datasets for each error distribution; every differentially
private release uses fresh privacy noise in each Monte Carlo replication.

Checkpoint files are created in `simulation_results/` during full runs. They
allow an interrupted run to resume. Delete only the relevant checkpoint
directory when intentionally starting the corresponding experiment from
scratch.

## Real-data application

The application uses the official 2024 ACS one-year California person-level
PUMS archive:

<https://www2.census.gov/programs-surveys/acs/data/pums/2024/1-Year/csv_pca.zip>

To download/cache the archive and construct the analytic cohort:

```text
Rscript --vanilla real_data_application.R prepare
```

To run the full analysis using 20 independent privacy-noise draws at each
privacy budget:

```text
Rscript --vanilla real_data_application.R all 20
```

The equivalent analysis command is:

```text
Rscript --vanilla real_data_application.R analyze 20
```

The `analyze` command automatically prepares the cohort when the cached cohort
does not exist. To regenerate figures from existing result files:

```text
Rscript --vanilla real_data_application.R figures
```

The downloaded archive and prepared cohort are cached in `data/`. Tables,
CSV files, RDS output, and figures are written to `results/`. The principal
figures are `acs_descriptive.pdf` and `acs_method_comparison.pdf`.

## RStudio

Clicking **Source** only defines the functions; it does not start a long run.
After sourcing, the equivalent console calls are:

```r
simulation_main(c("both", "full", "6", "200"))
real_data_main(c("all", "20"))
```

## Reproducibility and runtime

The random seeds, parameter grids, Gaussian kernel, bandwidth rules,
initialization schemes, privacy parameters, and output construction are set
inside the scripts. Runtime depends strongly on the processor and number of
worker processes. The 200-replication simulation is computationally intensive
and can require many hours; the five-replication mode should be used first to
check the local installation.

The ACS archive is approximately 278 MB. Internet access is needed only for the
first download; subsequent runs use the cached copy.

