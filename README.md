# WetDry-metrics

Compilation and analysis of wet lab and dry lab metrics from duplex-seq experiments performed in the BBGlab

## Create an environment with the required packages

```{bash}
git clone https://github.com/bbglab/wetdry-metrics.git
cd wetdry-metrics
conda env create -f environment.yml
conda activate metrics-env
```

## How to run it

```{bash}
cd scripts
python BuildWetDryMetrics.py --runs_list runs_list.json \
            --output_data_dir <output_directory> \
            --wetlab_qc_metrics_file wetlab_qc_metrics.xlsx \
            --plot_dir <path_to_plot_folder>
```

`--wetlab_qc_metrics_file` flag is optional depending on the availability of library prep. metrics.

`--plot_dir` flag automatically triggers the plotting of some summary metrics by batch and is also optional.

## Entire procedure description

1. (optional; highly recommended) Collect the library prep metrics in the standardized format provided via the template (wetlab_qc_metrics.xlsx) available in the GitHub repository. (https://github.com/bbglab/wetdry-metrics)

2. Clone the GitHub repository: bbglab/WetDryMetrics and move inside the directory.

```bash
git clone https://github.com/bbglab/wetdry-metrics.git 
cd wetdry-metrics
```

Update the `runs_list.json` file, also provided as a template, with the new information from your last run. Add a name of the run and the path to its’ deepUMIcaller output directory.

3. Create an environment with the required dependencies:

```bash
conda env create -f environment.yml
conda activate metrics-env
```

4. Run the metrics compilation script with the following command:

```bash
python BuildWetDryMetrics.py --runs_list runs_list.json \
--output_data_dir <output_directory> \
--wetlab_qc_metrics_file wetlab_qc_metrics.xlsx \
--plot_dir <path_to_plot_folder>
```

`--wetlab_qc_metrics_file` flag is optional depending on the availability of library prep. metrics.
`--plot_dir` flag automatically triggers the plotting of some summary metrics by batch and is also optional.

Check the printed progress for details on the steps executed.

5. Check the main metrics as described in table below to ensure that the run was successful. More info on anticipated results.

| Metric | Type of metric | Description |
|---|---|---|
| WetLab>>Input (ng) | Recovery | Initial DNA amount used for library prep in ng. Recoveries are calculated based on this initial input value. |
| WetLab>>DNA after ligation (ng) | Recovery | DNA amount recovered after ligation in ng. Used to calculate mass recovery at ligation step. |
| WetLab>>fmol to PCR1 | Recovery | Number of molecules in fmol measured by qPCR added to the PCR1 step; an estimate of the initial diversity of the library prior to PCR1 amplification. |
| WetLab>>Recovery Input to ligation | Recovery | Mass recovery of DNA from input to ligation. Considers total mass recovery including both ligated and non-ligated DNA. |
| WetLab>>Recovery Input to ligation qPCR | Recovery | Recovery of ligated DNA measured by qPCR. |
| DryLab>>On-target unique molecules percentage | Library prep. worked | Percentage of unique molecules that map to the capture panel regions. Typically expected >70%, depending on panel size and capture steps. |
| DryLab>>Raw to DSC | Library prep. worked | Total number of unique molecules divided by the number of double-strand consensus (DSC) families; approximates how many raw reads on average are required to form a DSC. |
| DryLab>>SSC to DSC | Library prep. worked | Number of single-strand consensus (SSC) families divided by the number of double-strand consensus (DSC) families. |
| DryLab>>Family size | Sequencing estimation | Peak of the distribution of family sizes (number of sequencing reads of the same unique molecule that form a family). |
| DryLab>>GBs analysed | Sequencing estimation | Gigabytes of sequencing data analysed. |
| DryLab>>Total GBs for optimal | Sequencing estimation | Estimated GBs required to reach optimal sequencing (defined as ~18x amplification per on-target molecule). |
| DryLab>>GBs (analysed - optimal) | Sequencing estimation | Difference between analysed GBs and GBs required for optimal sequencing (analysed - optimal). |
| DryLab>>Percentage of duplicates on-target | Summary stats | Percentage of reads that are duplicates of on-target unique molecules. |
| DryLab>>Amplification bias on/off | Sequencing estimation | Ratio of on-target molecules amplified versus off-target molecules; higher ratio indicates greater selectivity for on-target sequences. |
| DryLab>>Depth | Summary stats | Duplex sequencing depth in the provided target regions (note: boundary padding may cause this to be an underestimation). |
| Combined>>Unique molecules sequenced vs qPCR | Sequencing estimation | Ratio of unique molecules observed versus the qPCR estimate. Values close to 1 indicate good qPCR diversity estimation. |
| Combined>>Recovery input to depth | Recovery | Recovered duplex depth normalized by input DNA quantity; estimates overall yield and protocol efficiency. |
