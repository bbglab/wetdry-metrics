#!/usr/local/bin/python

import pandas as pd
import click

from BuildDryLabMetricsTable import compute_drylab_metrics_table
from PlotMetricsPerCohort import *
from utils_definitions import *

def parse_and_calculate_wetlab_recovery(qc_metrics):

    # 1. Mass-based recovery from input to ligation calculation
    qc_metrics["Recovery Input to ligation"] = qc_metrics["DNA after ligation (ng)"] / qc_metrics["Input (ng)"]

    # 2. qPCR-based recovery calculation
    # only if Tapestation peak size (after ligation) is provided
    if 'TapeStation (after ligation)' in qc_metrics.columns :
        # Calculate fragment size excluding size of the adapter
        qc_metrics["Ligation nonadapter size (bp)"] = qc_metrics["TapeStation (after ligation)"] - adapter_size

        # Convert fmol to molecule count
        qc_metrics["qPCR unique molecules"] = qc_metrics['fmol to PCR1'] * avogadro_in_f_scale

        # Estimate DNA mass from unique molecules and size
        qc_metrics["Lig_qpcr_DNA"] = qc_metrics["qPCR unique molecules"] * \
                                                    qc_metrics["Ligation nonadapter size (bp)"] * bp_weight_in_ng
        
        # Calculate recovery                             
        qc_metrics["Recovery Input to ligation qPCR"] = qc_metrics["Lig_qpcr_DNA"] / qc_metrics["Input (ng)"]

    # 3. Column ordering and prefixing
    id_cols =  ['BATCH', 'WetLab ID', 'DryLab ID']
    qc_metrics = qc_metrics[ id_cols + [x for x in qc_metrics.columns if x not in id_cols]]
    qc_metrics.columns = id_cols + [x for x in qc_metrics.columns[3:] ]
    qc_metrics.columns = id_cols + [f'WetLab>>{x}' for x in qc_metrics.columns[3:] ]

    return qc_metrics

@click.command()
@click.option(
    '--runs_list',
    required=True,
    type=click.Path(exists=True),
    help='Path to the JSON file containing the list of runs.'
)
@click.option(
    '--wetlab_qc_metrics_file',
    required=False,
    help='OPTIONAL: Path to the wet lab QC metrics Excel file.'
)
@click.option(
    '--output_data_dir',
    required=True,
    help='Directory where processed data files will be stored.'
)
@click.option(
    '--plot_dir',
    required=False,
    help='OPTIONAL: Directory to save generated plots.'
)

def build_consolidated_metrics_table(runs_list, 
                                     wetlab_qc_metrics_file, 
                                     output_data_dir,
                                     plot_dir):
    
    # 1. Generate metrics from deepUMI runs (dry metrics)
    drylab_metrics = compute_drylab_metrics_table(runs_list,
                                                  output_data_dir)
    
    
    if wetlab_qc_metrics_file is None:

        click.echo(f'\nNot wetlab metrics file provided, only Drylab metrics table were generated and saved in {output_data_dir}')

        # Save drylab metrics table
        drylab_metrics = drylab_metrics[informative_columns_drylab]
        drylab_metrics.to_csv(f'{output_data_dir}/DryLabMetrics.tsv', sep = '\t', header = True, index = False)

        # OPTIONAL: if plot_dir provided, plot metrics in pdf
        if plot_dir :
            plot_metrics(drylab_metrics, 
                        output_plot_dir=plot_dir)
            
            click.echo(f'\n\nMetrics plots saved in {plot_dir}')

        return
    
    else :

        click.echo('\nProcessing wetlab metrics file to generate consolidated wetdry metrics table')

        # 2. Load and parse wetlab QC metrics table and calculate recoveries
        qc_wetlab_metrics = pd.read_excel(wetlab_qc_metrics_file, skiprows=0)
        parsed_wetlab_metrics = parse_and_calculate_wetlab_recovery(qc_wetlab_metrics)

        # 3. Merge wetlab and drylab metrics
        consolidated_metrics = parsed_wetlab_metrics.merge(drylab_metrics,
                                                    left_on='DryLab ID',
                                                    right_on='SAMPLE_ID',
                                                    how='outer')

        # 4. Compute an estiamtion of the number of unique molecules on target measured by qPCR and compute ratio of unique molecules sequenced vs expected by qPCR
        consolidated_metrics['WetLab>>qPCR unique molecules'] = consolidated_metrics['WetLab>>fmol to PCR1'] * \
                                                                    consolidated_metrics['DryLab>>PANEL2GENOMEratio'] * avogadro_in_f_scale 
        consolidated_metrics['Combined>>Unique molecules sequenced vs qPCR'] = consolidated_metrics['FamMetrics>>on_target.unique_molecules'] / consolidated_metrics['WetLab>>qPCR unique molecules']

        # 5. Calculate Recovery input to duplex depth
        consolidated_metrics[f'Combined>>Recovery input to depth'] = consolidated_metrics[f"DryLab>>Depth"] / \
                                                                            (consolidated_metrics['WetLab>>Input (ng)'] / genome_weight)

        # 6. Select informative columns and store table
        consolidated_metrics_informative = consolidated_metrics[informative_columns]
        consolidated_metrics_informative.to_csv(f'{output_data_dir}/WetDryMetrics.tsv', sep = '\t', header = True, index = False)

        click.echo(f'\n\nConsolidated WetDry metrics table done and stored in {output_data_dir}')

        # 7. OPTIONAL: if plot_dir provided, plot metrics in pdf

        if plot_dir :
            plot_metrics(consolidated_metrics, 
                        output_plot_dir=plot_dir)
            
            click.echo(f'\n\nMetrics plots saved in {plot_dir}')


if __name__ == "__main__":
    build_consolidated_metrics_table()
