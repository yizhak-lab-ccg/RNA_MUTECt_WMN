from collect_features import *
from train import *
from preprocess import *
from create_features import *


def create_somatics_maf(somatics):
    # create vcf from somatics
    tmp = somatics[
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
    # create the .vcf, add vcf headers, and save the vars in it
    vcf_to_maf = os.path.abspath(cancer_dir + '/somatics.vcf')
    with open(vcf_to_maf, 'w') as file:
        file.write('##fileformat=VCFv4.0\n')
        for i in range(1, 23):
            file.write(f"##contig=<ID={i}>\n")
        file.write("##contig=<ID=X>\n")
        file.write("##contig=<ID=Y>\n")
    vcf.to_csv(vcf_to_maf, mode='a', index=False, sep='\t')
    # run funcotator to get MAF file
    funco_datasource = os.path.abspath(f"{funcotator_dir}funcotator_dataSources.v1.7.20200521g")
    funco_output = os.path.abspath(f"{results_path}/somatics.maf")
    # run funcotator
    os.system(f"{tools}gatk Funcotator --variant {vcf_to_maf} --reference {reference_hg19} --ref-version hg19 "
              f"--data-sources-path {funco_datasource} --output {funco_output} --output-file-format MAF "
              f"--java-options '-Xmx6G' --force-b37-to-hg19-reference-contig-conversion --QUIET true")


def test_models(test_set):
    os.system(f"rm -f {results_path}/test*")  # clear old results
    X = test_set[features]
    y = test_set['is_real_keep']
    sum_pred = np.zeros(len(y))
    for model_file in os.listdir(models_path):
        model = pickle.load(open(os.path.abspath(models_path + model_file), 'rb'))
        pred = model.predict(X)
        sum_pred += pred

    test_set['pred'] = np.where(sum_pred > num_folds / 2, 1, 0)  # choose label by majority vote

    # creating a MAF file (Funcotator output) from variants classified as somatic
    somatics = test_set[test_set['pred'] == 1]
    create_somatics_maf(somatics)


def build_and_train():
    preprocess()
    create_features()
    collect_features()
    train_set, test_set = prepare_train_test()
    print("Data is ready for training")
    train_models(train_set)
    print("finished training")
    return test_set


def main():
    test_set = build_and_train()
    test_models(test_set)


if __name__ == "__main__":
    main()
