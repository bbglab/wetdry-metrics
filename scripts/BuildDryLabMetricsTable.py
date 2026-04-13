#!/usr/local/bin/python

import pandas as pd
import json

from utils_definitions import *

#### Define VARIABLES ####

quality_to_shortname = {'raw_qc'    : 'raw',
                        'allm_qc'   : 'am',
                        'duplex_qc' : 'duplex'
                        }
#### FUNCTIONS ####

# FastQC 

def load_fastqc_data(fastqc_general_stats, 
                     columns_of_interest = [
                                            'percent_gc', 'avg_sequence_length', 'median_sequence_length',
                                            'total_sequences', 'percent_duplicates', 'percent_fails'
                                            ],
                     parallel = False
                     ):
    fastqc_data = pd.DataFrame(fastqc_general_stats).T.reset_index()
    fastqc_data[["sample", "read_pair"]] = [('_'.join(x.split('_')[:-1]), x.split('_')[-1]) for x in fastqc_data["index"]]

    # TODO here we could improve the handling of the information
    # we have removed the handling of individual read information
    output_data = fastqc_data.groupby("sample").agg({
        'percent_gc': 'mean',
        'avg_sequence_length': 'mean',
        'median_sequence_length': 'mean',
        'total_sequences': 'sum',
        'percent_duplicates': 'mean',
        'percent_fails': 'mean'
    }).reset_index()
    output_data = output_data[["sample"] + columns_of_interest]
    output_data.columns = ["sample"] + [f'FastQC>>{x}' for x in columns_of_interest]

    # not sure this is needed
    output_data["sample"] = output_data["sample"].astype(str)

    if parallel:
        output_data = fastqc_from_parallel(output_data)

    return output_data


def fastqc_from_parallel(fastqc_info):
    """
    When running par_chr, fastqc is run on each part, so we need to aggregate the data
    """
    # split sample column in two columns 'sample_id' and part
    fastqc_info['part'] = [sid.split('_')[-1] for sid in fastqc_info['sample']]
    fastqc_info['new_sampleid'] = [
        '_'.join(sid.split('_')[:-1]) if 'LPART' in part else sid
        for sid, part in zip(fastqc_info['sample'], fastqc_info['part'])
    ]

    fastqc_info['duplicates'] = fastqc_info['FastQC>>total_sequences'] * fastqc_info['FastQC>>percent_duplicates'] / 100

    # FIXME: percent_gc would need to be recalculated as a weighted mean, but for now we will take the median
    fastqc_info = fastqc_info[[
        'new_sampleid', 'FastQC>>total_sequences',
        'FastQC>>percent_gc', 'part', 'duplicates'
    ]]

    fastqc_info = fastqc_info.groupby('new_sampleid').agg({
        'FastQC>>total_sequences': 'sum',
        'FastQC>>percent_gc': 'median',
        'duplicates': 'sum'
    }).reset_index()

    fastqc_info['FastQC>>percent_duplicates'] = round(
        fastqc_info['duplicates'] / fastqc_info['FastQC>>total_sequences'] * 100, 3
    )

    fastqc_info['sample'] = fastqc_info['new_sampleid']
    fastqc_info = fastqc_info[['sample', 'FastQC>>total_sequences', 'FastQC>>percent_duplicates', 'FastQC>>percent_gc']]

    return fastqc_info


# Qualimap

def get_qualimap_data(
    qualimap_general_stats,
    quality,
    columns_of_interest=['total_reads', 'mapped_reads',
                         'general_error_rate', 'mean_coverage', 'regions_size',
                         'regions_mapped_reads', 'percentage_aligned',
                         'percentage_aligned_on_target', 'median_coverage',
                         'median_insert_size']
):
    qualimap_all = pd.DataFrame()

    for sample in qualimap_general_stats.keys():
        qualimap_samp = pd.DataFrame.from_dict(qualimap_general_stats[sample], orient='index').T
        qualimap_samp["sample"] = sample
        qualimap_samp = qualimap_samp[['sample'] + columns_of_interest]

        qualimap_samp.columns = ['sample'] + [f'{quality}.{x}' for x in columns_of_interest]

        qualimap_all = pd.concat((qualimap_all, qualimap_samp))
        qualimap_all["sample"] = qualimap_all["sample"].astype(str)

    return qualimap_all


# Family size metrics

def fam_metrics2df(file, prefix='all'):
    data = pd.read_csv(file, sep='\t')

    # Convert relevant columns to numeric, using errors='coerce' to handle non-numeric values
    numeric_columns = ['raw_reads', 'duplicates', 'sscs', 'dscs',
                       'expected_dscs', 'recovery_of_dscs', 'unique_reads', 'uq_reads_duplex',
                       'unique_molecules', 'raw_x_dscs', 'raw_x_sscs', 'sscs_x_dscs',
                       'duplex_raw_reads', 'duplex_sscs', 'duplex_raw_x_dscs',
                       'duplex_raw_x_sscs', 'noduplex_raw_reads', 'noduplex_sscs',
                       'noduplex_raw_x_sscs', 'peak_size']

    # Convert the specified columns to numeric
    for col in numeric_columns:
        if col in data.columns:  # Check if column exists before conversion
            data[col] = pd.to_numeric(data[col], errors='coerce')

    # transform multiple rows per sample to single row per sample
    pivot_df = data.pivot_table(index=['sample', 'raw_reads', 'sscs', 'raw_x_sscs', 'duplicates'], columns='quality', aggfunc='first')
    pivot_df.columns = [f'{col[1]}.{col[0]}' for col in pivot_df.columns]

    for prune in ["unique_molecules", "unique_reads"]:
        pivot_df[prune] = pivot_df[f"allm.{prune}"].copy()
        pivot_df = pivot_df.drop([f"allm.{prune}", f"duplex.{prune}"], axis='columns')

    data = pivot_df.reset_index()
    data.columns = ["sample"] + [f'FamMetrics>>{prefix}.{x}' for x in data.columns[1:]]
    data["sample"] = data["sample"].astype(str)

    return data


# Family size metrics - off target metrics

def compute_off_targets_metrics(data_family):
    """
    With the information from on target and total reads compute the familymetrics for off targets.
    TODO: this could be included in deepUMIcaller, but for now it is easier to do it here.
    """

    off_target_df = data_family[["sample"]].copy()

    for column in ['sscs', 'raw_reads', 'unique_molecules', 'unique_reads']:
        off_target_df[f'FamMetrics>>off_target.{column}'] = (
            data_family[f'FamMetrics>>all.{column}'] - data_family[f'FamMetrics>>on_target.{column}']
        )

    off_target_df['FamMetrics>>off_target.raw_x_sscs'] = (
        off_target_df['FamMetrics>>off_target.raw_reads'] / off_target_df['FamMetrics>>off_target.sscs']
    )

    column = 'dscs'

    for level in ["allm", "duplex"]:
        off_target_df[f'FamMetrics>>off_target.{level}.{column}'] = (
            data_family[f'FamMetrics>>all.{level}.{column}'] - data_family[f'FamMetrics>>on_target.{level}.{column}']
        )

        off_target_df[f'FamMetrics>>off_target.{level}.raw_x_dscs'] = (
            off_target_df[f'FamMetrics>>off_target.raw_reads'] / off_target_df[f'FamMetrics>>off_target.{level}.dscs']
        )
        off_target_df[f'FamMetrics>>off_target.{level}.sscs_x_dscs'] = (
            off_target_df[f'FamMetrics>>off_target.sscs'] / off_target_df[f'FamMetrics>>off_target.{level}.dscs']
        )

    return off_target_df


# Compute FOR OPTIMAL metrics

def compute_desired_sequencing_metrics(family_data, desired_amplification=target_amplification):
    """
    Here we are computing the metrics that are used for estimating the optimal sequencing.

    TODO: we are not computing the distance from the GBs analyzed to the GBs optimal,
            it would be great to do it in this script, or even this function
    """

    family_data['DryLab>>On-target unique molecules percentage'] = (
        family_data['FamMetrics>>on_target.unique_molecules'] / family_data['FamMetrics>>all.unique_molecules']
    )
    family_data['DryLab>>Off-target unique molecules percentage'] = (
        family_data['FamMetrics>>off_target.unique_molecules'] / family_data['FamMetrics>>all.unique_molecules']
    )

    family_data['DryLab>>Amplification bias on/off'] = (
        (family_data['FamMetrics>>on_target.raw_reads'] / family_data['FamMetrics>>on_target.unique_molecules'])
        / (family_data['FamMetrics>>off_target.raw_reads'] / family_data['FamMetrics>>off_target.unique_molecules'])
    )

    # UM * %on * desired_amplification = reads for optimal on target
    # UM * %off * (desired_amplification / amplification_bias.on/off )  = reads for optimal off target

    family_data['FamMetrics>>molecules_for_optimal_on_target'] = (
        family_data['FamMetrics>>on_target.unique_molecules'] * desired_amplification
    )

    family_data['FamMetrics>>molecules_for_optimal_off_target'] = (
        family_data['FamMetrics>>off_target.unique_molecules'] * (desired_amplification / family_data['DryLab>>Amplification bias on/off'])
    )

    family_data['FamMetrics>>total_molecules_for_optimal'] = (
        family_data['FamMetrics>>molecules_for_optimal_on_target'] + family_data['FamMetrics>>molecules_for_optimal_off_target']
    )
    family_data['DryLab>>Total GBs for optimal'] = (
        family_data['FamMetrics>>total_molecules_for_optimal'] * 300 / 1000000000
    )

    return family_data

# Compute metrics from family metrics

def compute_metrics_from_familymetrics(path_family_metrics):
    data = fam_metrics2df(f"{path_family_metrics}/metrics_summary.tsv")
    data_ont = fam_metrics2df(f"{path_family_metrics}ontarget/metrics_summary.tsv", prefix='on_target')

    data_family = data.merge(data_ont, on='sample')

    off_targets_families = compute_off_targets_metrics(data_family)
    data_family_joined = data_family.merge(off_targets_families, on='sample')

    return compute_desired_sequencing_metrics(data_family_joined)


def compute_drylab_metrics_table(runs_list, output_data_dir):

    #### Load runs from JSON file
    with open(runs_list, "r") as file:
        data = json.load(file)

    # Convert JSON data to the exact structure: list of tuples
    runs_list_missing = [(entry["name"], entry["path"], entry.get('parallel', False)) for entry in data]

    # Print the result for validation
    for run in runs_list_missing:
        print(f"Name: {run[0]}, Path: {run[1]}, Parallel : {run[2]}")

    # Specify the metrics path
    metrics_relative_path_options = ['metrics/duplex/familymetrics', 'familymetrics']

    # Generate metrics for missing runs
    all_data = pd.DataFrame()

    for i, path_run, split_fq in runs_list_missing:
        list_sample_qualities = ['raw_qc', 'allm_qc']
        # print(i)
        data_family_complete = None
        last_exception = None
        for metrics_relative_path in metrics_relative_path_options:
            try:
                path_metrics = f"{path_run}/{metrics_relative_path}"
                data_family_complete = compute_metrics_from_familymetrics(path_metrics)
                # Successfully loaded metrics; no need to try other paths
                break
            except Exception as e:
                # Record the exception and try the next possible metrics path
                last_exception = e
                continue

        if data_family_complete is None:
            # All metrics paths failed for this run; report and skip to next run
            print(
                f"Error occurred while obtaining the family metrics for {i}.\n"
                f"Tried paths: {metrics_relative_path_options}.\n"
                f"Last error: {last_exception}"
            )
            continue

        #### Define info columns : TISSUE, ANALYSIS BATCH, PATH
        data_family_complete["Run ID"] = '_'.join(i.split("_")[:2])
        data_family_complete["PATH"] = path_run

        #### Rename columns 
        data_family_complete = data_family_complete.rename(columns={
            'FamMetrics>>on_target.duplex.raw_x_dscs': 'DryLab>>Raw to DSC',
            'FamMetrics>>on_target.duplex.sscs_x_dscs': 'DryLab>>SSC to DSC',
            'FamMetrics>>on_target.duplex.peak_size' : 'DryLab>>Family size',
            'FamMetrics>>on_target.duplicates' : 'DryLab>>Percentage of duplicates on-target',
        })        

        #### Load multiqc results
        path_multiqc = f"{path_run}/multiqc/multiqc_data"

        with open(f"{path_multiqc}/multiqc_data.json", 'r') as file:
            dict_data = json.load(file)

        # Parse dictionary with multiqc results to retrieve general stats
        header_2_index = {}

        keys_multiqc_dict = list(dict_data['report_general_stats_headers'].keys())
        keys_names = {"qualimap_raw": "raw_qc",
                     "qualimap_all": "allm_qc",
                     "qualimap_duplex": "duplex_qc",
                     "fastqc": "fastqc"}

        for key in keys_multiqc_dict:
            header_2_index[keys_names[key]] = key

        general_stats = dict_data['report_general_stats_data']

        #### Load bamqc available raw and allm data safely
        list_sample_qualities = ['raw_qc', 'allm_qc']
        bamqc_raw = get_qualimap_data(general_stats[header_2_index['raw_qc']], quality_to_shortname['raw_qc'])

        #### if parallel option was used, bamqc data needs to be compile - only for raw
        # Sanity checks : parallel option
        if not par_chr and bamqc_raw['sample'].str.contains('LPART').any():
            raise ValueError(
                f"In {i}, parallel is set to False but sample names indicate parallel chromosome processing. "
                "Please check the input data."
            )

        if par_chr and not bamqc_raw['sample'].str.contains('LPART').any():
            raise ValueError(
                f"In {i}, parallel is set to True but sample names indicate parallel chromosome processing WAS NOT applied. "
                "Please check the input data."
            )

        if par_chr:

            # split sample column in two columns 'sample_id' and part
            bamqc_raw['part'] = [sid.split('_')[-1] for sid in bamqc_raw['sample']]
            bamqc_raw['new_sampleid'] = [
                '_'.join(sid.split('_')[:-1]) if 'LPART' in part else sid
                for sid, part in zip(bamqc_raw['sample'], bamqc_raw['part'])
            ]

            bamqc_raw = bamqc_raw[[
                'new_sampleid', 'raw.total_reads', 'raw.mapped_reads',
                'raw.general_error_rate', 'raw.mean_coverage', 'raw.regions_size',
                'raw.regions_mapped_reads', 'raw.median_coverage',
                'raw.median_insert_size', 'part'
            ]]

            # group by sample_id and quality and sum the columns that can be summed
            bamqc_raw = bamqc_raw.groupby(['new_sampleid']).agg({
                'raw.total_reads': 'sum',
                'raw.mapped_reads': 'sum',
                'raw.general_error_rate': 'median',
                'raw.mean_coverage': 'median',
                'raw.regions_size': 'first',
                'raw.regions_mapped_reads': 'sum',
                'raw.median_coverage': 'median',
                'raw.median_insert_size': 'median'
            }).reset_index()

            # calculate raw.percentage_aligned and raw.percentage_aligned_on_target
            bamqc_raw['raw.percentage_aligned'] = bamqc_raw['raw.mapped_reads'] / bamqc_raw['raw.total_reads'] * 100
            bamqc_raw['raw.percentage_aligned_on_target'] = bamqc_raw['raw.regions_mapped_reads'] / bamqc_raw['raw.total_reads'] * 100

            bamqc_raw['sample'] = bamqc_raw['new_sampleid']
            bamqc_raw = bamqc_raw[[
                'sample', 'raw.total_reads', 'raw.mapped_reads',
                'raw.general_error_rate', 'raw.mean_coverage', 'raw.regions_size',
                'raw.regions_mapped_reads', 'raw.percentage_aligned',
                'raw.percentage_aligned_on_target', 'raw.median_coverage',
                'raw.median_insert_size'
            ]]

        bamqc_allm = get_qualimap_data(general_stats[header_2_index['allm_qc']], quality_to_shortname['allm_qc'])

        bamqc_duplex = get_qualimap_data(general_stats[header_2_index['duplex_qc']], quality_to_shortname['duplex_qc'])
        list_sample_qualities.append('duplex_qc')

        # Merge available data
        bam_qc_data_ra = bamqc_raw.merge(bamqc_allm, on='sample', how='left')
        bam_qc_data_ra = bam_qc_data_ra.merge(bamqc_duplex, on='sample', how='left')

        bam_qc_data_ra.columns = ['sample'] + [f'BamQC>>{x}' for x in bam_qc_data_ra.columns[1:]]

        # Rename depth column
        bam_qc_data_ra = bam_qc_data_ra.rename(columns={'BamQC>>duplex.mean_coverage': 'DryLab>>Depth'})

        #### Load FastQC data
        fastqc_info = load_fastqc_data(
            general_stats[header_2_index['fastqc']],
            ['total_sequences', 'percent_duplicates', 'percent_gc'],
            parallel=par_chr
        )

        # Translate analysed data to Gbs
        fastqc_info['DryLab>>GBs analysed'] = fastqc_info['FastQC>>total_sequences'] * read_size / 1e9

        data_fm_bamqc = data_family_complete.merge(bam_qc_data_ra, on='sample', how='left')
        data_fm_bamqc_fq = data_fm_bamqc.merge(fastqc_info, on='sample', how='left')

        all_data = pd.concat((all_data, data_fm_bamqc_fq))

    if all_data.shape[0] > 0:
        all_data = all_data.reset_index(drop=True)
        # change column names
        all_data = all_data[["sample", "Run ID", "PATH"] + [x for x in all_data.columns if x not in ["sample", "Run ID", "PATH"]]]
        all_data.columns = ["SAMPLE_ID"] + list(all_data.columns[1:])

        #### Compute sequence GBs missing to optimal
        all_data['DryLab>>GBs (analyzed-optimal)'] = all_data['DryLab>>GBs analysed'] - all_data['DryLab>>Total GBs for optimal']

        #### Add columns with panel info and panel2genome ratio
        all_data["DryLab>>PANEL2GENOMEratio"] = all_data["BamQC>>raw.regions_size"].copy() / size_whole_genome

        #### Calculate GBs analyzed to optimal ratio
        all_data['DryLab>>GBs_analysed/optimal'] = all_data['DryLab>>GBs analysed'] / all_data["DryLab>>Total GBs for optimal"]

    print(f'\n\nDryLab metrics table successfully saved to : {output_data_dir}')

    return all_data



