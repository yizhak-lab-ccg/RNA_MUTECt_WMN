import pandas as pd
import numpy as np
import os
from os import path
from pathlib import Path
from capy import mut
from config import *


def extract_af(all_vars):
    print("extracting af from bcftools output...")
    for file_name in os.listdir(cancer_dir + "/BCF_TOOLS/final_outputs"):
        if "parsed" in file_name:
            continue
        cur_output_file_path = os.path.abspath(cancer_dir + "/BCF_TOOLS/final_outputs/" + file_name)
        db_type = os.path.splitext(file_name)[0].split("_")[0]  # assumption: file name = <db_type>_0002.txt
        col_names = ['#CHROM', 'POS', 'ALT', 'INFO']
        with open(cur_output_file_path, "r") as output_file:
            lines = output_file.readlines()
        if len(lines) == 0:  # output file is empty, move to the next file
            continue
        for index, line in enumerate(lines):
            if line.strip("\n").startswith("#CHROM"):
                break  # now index = row number of the real header
        df = pd.read_table(cur_output_file_path, usecols=col_names, dtype={'#CHROM': 'str'}, skiprows=index)
        df = df.rename(columns={'#CHROM': 'Chromosome', 'POS': 'position', 'ALT': 'alt_allele'})
        AF_symbol = 'AF'  # if db_type in ['gnomad', 'genome1000']:
        if db_type == 'dbsnp':
            AF_symbol = 'CAF'
        if db_type == 'esp':
            AF_symbol = 'MAF'
        info_params = (
            df["INFO"].str.split(";", expand=True).stack().str.split("=", expand=True).reset_index(level=1, drop=True))
        add_params = pd.concat([df, info_params.groupby([info_params.index, info_params[0]])[1].sum().unstack()],
                               axis=1).drop("INFO", axis=1)
        add_params = add_params[['Chromosome', 'position', 'alt_allele', AF_symbol]]
        if db_type == "dbsnp":
            add_params[AF_symbol] = add_params[AF_symbol].str.split(',').str[1]
        elif db_type == 'esp':
            add_params[AF_symbol] = add_params[AF_symbol].str.split(',').str[2]
        af_col = f'{db_type}_af'
        add_params = add_params.rename(columns={AF_symbol: af_col})
        all_vars = pd.merge(all_vars, add_params, how='left')
        all_vars[db_type] = np.where(all_vars[af_col].isna(), 0, 1)  # binary db column
        # replacing nan af values with means
        all_vars[af_col] = pd.to_numeric(all_vars[af_col])
        af_mean = np.nanmean(all_vars[af_col], axis=0)
        all_vars[af_col] = np.where(all_vars[af_col].isna(), af_mean, all_vars[af_col])
        print(f"finished extracting af from {db_type}")
    return all_vars


def calc_pon(capy_input, sample_type, binary, ref, all_vars, output_path):
    capy_input[f"log_like_{sample_type}"] = mut.filter_mutations_against_token_PoN(M=capy_input, ponfile=binary,
                                                                                   ref=ref)
    pons_input = capy_input[['chr', 'pos']].drop_duplicates()
    pons = pd.DataFrame(mut.get_pon(M=pons_input, ponfile=binary, ref=ref))
    assert (pons.shape[0] == pons_input.shape[0])
    assert (pons.shape[1] == 8)
    for index in range(8):
        pons_input[f'pon_{sample_type}_{index + 1}'] = pons[index]
    full_pon_output = pd.merge(capy_input, pons_input)
    assert (full_pon_output.shape[0] == capy_input.shape[0])

    full_pon_output = full_pon_output.astype({"chr": str})
    full_pon_output['chr'] = np.where(full_pon_output['chr'] == '23', 'X', full_pon_output['chr'])
    full_pon_output['chr'] = np.where(full_pon_output['chr'] == '24', 'Y', full_pon_output['chr'])
    full_pon_output = full_pon_output.rename(
        columns={"chr": "Chromosome", "pos": "position", "n_ref": "n_ref_count", "n_alt": "n_alt_count"})

    chr_vars = all_vars[all_vars['Chromosome'] == '1']
    chr_pons = full_pon_output[full_pon_output['Chromosome'] == '1']
    chr_vars_with_pon = pd.merge(chr_vars, chr_pons, how='left')
    chr_vars_with_pon.to_csv(output_path, index=False)

    chrs = [str(x) for x in range(2, 23)]
    chrs += ['X', 'Y']
    for cur_chr in chrs:
        chr_vars = all_vars[all_vars['Chromosome'] == cur_chr]
        rows_count = chr_vars.shape[0]
        chr_pons = full_pon_output[full_pon_output['Chromosome'] == cur_chr]
        chr_vars_with_pon = pd.merge(chr_vars, chr_pons, how='left')
        new_rows_count = chr_vars_with_pon.shape[0]
        if new_rows_count != rows_count:
            print(f"wrong rows in chr {cur_chr}!")
        chr_vars_with_pon.to_csv(output_path, mode='a', index=False, header=False)


def extract_pon(all_vars):
    print("calculating pon")
    Path(pon_dir).mkdir(parents=True, exist_ok=True)
    # need to rename columns to use CApy methods
    capy_input = all_vars[['Chromosome', 'position', 'n_ref_count', 'n_alt_count']]
    capy_input = capy_input.drop_duplicates()
    capy_input = capy_input.rename(
        columns={"n_ref_count": "n_ref", "n_alt_count": "n_alt", "Chromosome": "chr", "position": "pos"})

    # converting chromosomes names to numbers
    capy_input['chr'] = np.where(capy_input['chr'] == 'X', '23', capy_input['chr'])
    capy_input['chr'] = np.where(capy_input['chr'] == 'Y', '24', capy_input['chr'])
    capy_input = capy_input.astype({"chr": int})
    all_vars_RNA_output_path = os.path.abspath(cancer_dir + "/all_vars_after_RNA_pon.csv")
    all_vars_DNA_output_path = os.path.abspath(cancer_dir + "/all_vars_after_DNA_pon.csv")
    calc_pon(capy_input, 'RNA', RNA_pon_binary, reference_hg19, all_vars, all_vars_RNA_output_path)
    all_vars_RNA_output = pd.read_csv(all_vars_RNA_output_path, dtype={'Chromosome': 'str'})
    calc_pon(capy_input.drop(columns=['log_like_RNA']), 'DNA', DNA_pon_binary, reference_hg19, all_vars_RNA_output,
             all_vars_DNA_output_path)


def extract_funco_data(all_vars):
    funco_output_path = os.path.abspath(cancer_dir + "/funco_output.vcf")
    with open(funco_output_path, "r", encoding='latin-1') as output_file:
        lines = output_file.readlines()
    for index,line in enumerate(lines):
        if line.strip("\n").startswith("##INFO"):
            funco_features = line.split("|")
        elif line.strip("\n").startswith("#CHROM"):
            break
    funco_features[0] = funco_features[0].split(": ")[1]
    funco_features[-1] = funco_features[-1].replace('">', '')
    funco_output = pd.read_table(funco_output_path, encoding='latin-1', dtype={'#CHROM': 'str'}, skiprows=index)
    funco_output = funco_output.rename(columns={'#CHROM': 'Chromosome', 'POS': 'position', 'ALT': 'alt_allele'})
    funco_output["Chromosome"] = funco_output["Chromosome"].str.replace("chr", "")
    funco_output["INFO"] = funco_output["INFO"].str.replace("FUNCOTATION=\[", "")
    funco_output["INFO"] = funco_output["INFO"].str.replace("\]", "")
    funco_output[funco_features] = funco_output['INFO'].str.split('|', expand=True)
    funco_output = funco_output[
        ["Chromosome", "position", "alt_allele", 'Gencode_34_hugoSymbol', 'Gencode_34_variantClassification',
         'Gencode_34_variantType', 'Gencode_34_referenceContext', 'Gencode_34_gcContent']]
    funco_output = funco_output.rename(columns={'Gencode_34_hugoSymbol': 'Hugo_Symbol',
                                                'Gencode_34_variantClassification': 'Variant_Classification',
                                                'Gencode_34_variantType': 'Variant_Type',
                                                "Gencode_34_referenceContext": "ref_context",
                                                "Gencode_34_gcContent": "gc_content"})
    all_vars = pd.merge(all_vars, funco_output,
                        how='left')
    classification_to_remove = ['IGR', 'Intron', 'RNA', 'lincRNA']
    all_vars["classification_to_remove"] = np.where(all_vars['Variant_Classification'].isin(classification_to_remove),
                                                    1, 0)
    return all_vars


def collect_features():
    print("running feature collection...")

    input_all_vars_path = os.path.abspath(cancer_dir + '/all_vars_after_preprocess.csv')
    if not path.isfile(input_all_vars_path):
        print('Missing output of preprocess phase, please run again.')
        return
    all_vars = pd.read_csv(input_all_vars_path, dtype={'Chromosome': 'str'})

    all_vars_with_af = extract_af(all_vars)
    print("done af")
    all_vars_with_funco = extract_funco_data(all_vars_with_af)
    print("done funcotator")
    extract_pon(all_vars_with_funco)
    print("done pon")
