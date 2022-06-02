from config import *
from pathlib import Path
import os


def prepare_bcftools_input(bcftools_cancer_dir):
    raw_input_path = os.path.abspath(bcftools_cancer_dir + '/unique_variants.vcf')
    if not os.path.isfile(raw_input_path):
        if not os.path.isfile(f"{raw_input_path}.gz"):
            raise Exception("bcftools input is missing!")  # if there's also no vcf file, can't continue!
        else:
            return
    else:  # vcf exists, delete old versions and re-create them
        os.system(f"rm -f {raw_input_path}.gz")
        os.system(f"rm -f {raw_input_path}.gz.tbi")
        os.system(f"{tools}bgzip {raw_input_path}")  # create the .gz file
        os.system(f"{tools}tabix -p vcf {raw_input_path}.gz")  # create the index file


def download_bcftools_dbs():
    dbsnp_path = os.path.abspath(bcftools_dir + '/All_20180423.vcf.gz')
    dbsnp_link = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz"
    gnomad_path = os.path.abspath(bcftools_dir + '/gnomad.exomes.r2.1.1.sites.vcf.bgz')
    gnomad_link = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz"
    thousand_genome_path = os.path.abspath(bcftools_dir + '/ALL.2of4intersection.20100804.genotypes.vcf.gz')
    thousand_genome_link = "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz"
    esp_path = os.path.abspath(bcftools_dir + '/merged.vcf.gz')  # esp file is a required input
    dbs_files = {"dbsnp": [dbsnp_path, dbsnp_link], "gnomad": [gnomad_path, gnomad_link],
                 "thousand_genome": [thousand_genome_path, thousand_genome_link],
                 "esp": [esp_path]}
    for db, values in dbs_files.items():
        if db == "esp":
            continue
        if not os.path.isfile(values[0]):
            os.system(f"wget -P {bcftools_dir} {values[1]}")
        if not os.path.isfile(values[0] + ".tbi"):
            os.system(f"wget -P {bcftools_dir}  {values[1]}.tbi")
    print("all dbs are downloaded and ready")
    return dbs_files


def run_bcftools_on_db(db_path, db_name, bcftools_cancer_dir):
    print(f"starting bcftools on {db_name}. This will take a while...")
    output_dir = os.path.abspath(bcftools_cancer_dir + "/bcftools_output")
    input_path = os.path.abspath(bcftools_cancer_dir + '/unique_variants.vcf.gz')
    final_outputs = os.path.abspath(bcftools_cancer_dir + "/final_outputs")
    # use bcftools to intersect variants with current db
    os.system(f"{tools}bcftools isec -p {output_dir} -w1 -Oz {db_path} {input_path}")
    # unzip intersection output
    os.system(f"{tools}bgzip -d -c {output_dir}/0002.vcf.gz > {output_dir}/0002.vcf")
    os.system(f"mv {output_dir}/0002.vcf  {final_outputs}/{db_name}_0002.vcf")
    os.system(f"rm -f {output_dir}/*")


def run_bcftools():
    print("running bcftools...")
    # creating directories if needed
    bcftools_cancer_dir = os.path.abspath(cancer_dir + "/BCF_TOOLS/")
    Path(resource_dir + "/BCF_tools_dbs").mkdir(parents=True, exist_ok=True)
    Path(bcftools_cancer_dir + "/bcftools_output").mkdir(parents=True, exist_ok=True)
    Path(bcftools_cancer_dir + "/final_outputs").mkdir(parents=True, exist_ok=True)
    try:
        prepare_bcftools_input(bcftools_cancer_dir)
        print("bcftools input file is ready")
    except Exception as e:
        print(e)
        return
    dbs_files = download_bcftools_dbs()
    for db in dbs_files:
        run_bcftools_on_db(dbs_files[db][0], db, bcftools_cancer_dir)


def run_funcotator():
    print("running Funcotator...")
    # preparation for running
    Path(resource_dir + "/funcotator").mkdir(parents=True, exist_ok=True)
    if not os.path.isfile(f"{funcotator_dir}dataSources.v1.7.20200521g.tar.gz"):
        os.system(f"{tools}gatk FuncotatorDataSourceDownloader --germline --validate-integrity "
                  f"--extract-after-download --output {funcotator_dir}dataSources.v1.7.20200521g.tar.gz")
    vcf_path = os.path.abspath(f"{cancer_dir}unique_variants.vcf")
    funco_datasource = os.path.abspath(f"{funcotator_dir}funcotator_dataSources.v1.7.20200521g")
    funco_output = os.path.abspath(f"{cancer_dir}funco_output.vcf")
    # run funcotator
    os.system(f"{tools}gatk Funcotator --variant {vcf_path} --reference {reference_hg19} --ref-version hg19 "
              f"--data-sources-path {funco_datasource} --output {funco_output} --output-file-format VCF "
              f"--java-options '-Xmx6G' --force-b37-to-hg19-reference-contig-conversion --QUIET true")


def create_features():
    run_bcftools()
    run_funcotator()
