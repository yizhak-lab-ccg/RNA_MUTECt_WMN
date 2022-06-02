import os
import pandas as pd
import numpy as np
from config import *
from pathlib import Path


def create_df(path, file_type):
    assert (file_type in ["call_stats", "maf"])
    if file_type == "call_stats":
        col_names = ['contig', 'position', 'ref_allele', 'alt_allele', 'tumor_name',
                     'total_reads', 'tumor_power', 'contaminant_lod', 't_ref_count', 't_alt_count',
                     't_ref_sum', 't_alt_sum', 'n_ref_count', 'n_alt_count', 'n_ref_sum', 'n_alt_sum',
                     't_lod_fstar', 'tumor_f', 'judgement', 'failure_reasons', 'context']
        # read relevant columns from call_stats file
        df = pd.read_table(path, usecols=col_names, header=1, dtype={'contig': 'str'})
        df = df.rename(columns={"tumor_name": "Tumor_Sample_Barcode", "contig": "Chromosome"})

    if file_type == "maf":
        col_names = ["Chromosome", "Start_position", "Tumor_Sample_Barcode"]
        df = pd.read_table(path, usecols=col_names, dtype={'Chromosome': 'str'})
        df = df.rename(columns={"Start_position": "position"})

    return df


def filter_call_stats_df(df):
    df.failure_reasons = df.failure_reasons.fillna('')
    # remove results of GL chromosome
    df = df[~df.Chromosome.str.contains("GL")]
    # remove results of MT chromosome
    df = df[df.Chromosome != 'MT']
    # keep only non-fail or valid fail reasons (=germline)
    non_germ_failure = ["clustered_read_position", "fstar_tumor_lod", "nearby_gap_events", "seen_in_panel_of_normals",
                        "poor_mapping_region_alternate_allele_mapq", "poor_mapping_region_mapq0",
                        "possible_contamination",
                        "strand_artifact", "triallelic_site"]
    df = df[(df.failure_reasons == "") | ~(df.failure_reasons.str.contains('|'.join(non_germ_failure)))]
    return df


def mark_noise(all_vars):
    # only vars that also exist in the MAF files are not noise.
    # after this function, the value in the "is_real_keep" columns indicates if the var is real or noise:
    # 1 is real and 0 is noise.
    dfs = []
    for maf_file in os.listdir(maf_files_dir):
        cur_maf_path = os.path.abspath(maf_files_dir + maf_file)
        df = create_df(cur_maf_path, file_type="maf")
        dfs.append(df)
    all_mafs = pd.concat(dfs, ignore_index=True)
    all_mafs = all_mafs.drop_duplicates()
    merged = pd.merge(all_vars, all_mafs, how='left', indicator='Exist')
    all_vars['is_real_keep'] = np.where(merged.Exist == 'both', 1, 0)
    return all_vars


def process_col_stats():
    dfs = []
    for file in os.listdir(call_stats_files_dir):
        cur_file_path = os.path.abspath(call_stats_files_dir + file)
        df = create_df(cur_file_path, file_type="call_stats")
        df = filter_call_stats_df(df)
        dfs.append(df)
    all_vars = pd.concat(dfs, ignore_index=True)
    all_vars = all_vars.drop_duplicates()
    # converting chromosomes names to ints for sorting
    all_vars['Chromosome'] = np.where(all_vars['Chromosome'] == 'X', '23', all_vars['Chromosome'])
    all_vars['Chromosome'] = np.where(all_vars['Chromosome'] == 'Y', '24', all_vars['Chromosome'])
    all_vars = all_vars.astype({"Chromosome": int})
    # sorting
    all_vars = all_vars.sort_values(by=['Chromosome', 'position'], ignore_index=True)

    # converting back to strings
    all_vars = all_vars.astype({"Chromosome": str})
    all_vars['Chromosome'] = np.where(all_vars['Chromosome'] == '23', 'X', all_vars['Chromosome'])
    all_vars['Chromosome'] = np.where(all_vars['Chromosome'] == '24', 'Y', all_vars['Chromosome'])
    all_vars = mark_noise(all_vars)
    all_vars.to_csv(os.path.abspath(cancer_dir + "/all_vars_after_preprocess.csv"), index=False)
    return all_vars


def create_vcf(all_vars):
    tmp = all_vars[
        ['Chromosome', 'position', 'ref_allele', 'alt_allele']]
    tmp = tmp.drop_duplicates(subset=['Chromosome', 'position', 'alt_allele'])
    vcf = pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    vcf['#CHROM'] = tmp['Chromosome']
    vcf['POS'] = tmp['position']
    vcf['ID'] = '.'
    vcf['REF'] = tmp['ref_allele']
    vcf['ALT'] = tmp['alt_allele']
    vcf['QUAL'] = '.'
    vcf['FILTER'] = 'PASS'
    vcf['INFO'] = '.'

    Path(cancer_dir + "/BCF_TOOLS").mkdir(parents=True, exist_ok=True)
    # create the unique_variants.vcf, add vcf headers, and save the vars in it
    bcf_tools_file_path = os.path.abspath(cancer_dir + '/BCF_TOOLS/unique_variants.vcf')
    with open(bcf_tools_file_path, 'w') as file:
        file.write('##fileformat=VCFv4.0\n')
        for i in range(1, 23):
            file.write(f"##contig=<ID={i}>\n")
        file.write("##contig=<ID=X>\n")
        file.write("##contig=<ID=Y>\n")
    vcf.to_csv(bcf_tools_file_path, mode='a', index=False, sep='\t')
    os.system(f"cp {bcf_tools_file_path} {cancer_dir}unique_variants.vcf")


def preprocess():
    print("Preprocessing the input...")
    all_vars = process_col_stats()
    create_vcf(all_vars)
