# dirs and files configuration #
data_dir = "./Data"
cancer_type = ""  # must configure!
cancer_dir = f"{data_dir}/{cancer_type}/"
models_path = cancer_dir + "/models/"
results_path = cancer_dir + "/results/"
input_path = cancer_dir + "/input/"
call_stats_files_dir = input_path + "call_stats/"
maf_files_dir = input_path + "maf/"
resource_dir = data_dir + "/resource/"
pon_dir = resource_dir + "pon/"
ref_dir = resource_dir + "reference/"
bcftools_dir = resource_dir + "BCF_tools_dbs/"
funcotator_dir = resource_dir + "funcotator/"
# pon and ref file names must be configured! #
RNA_pon_binary = pon_dir + ""
DNA_pon_binary = pon_dir + ""
reference_hg19 = ref_dir + ""

assert(cancer_type != "")
assert (RNA_pon_binary != pon_dir and DNA_pon_binary != pon_dir)
assert (reference_hg19 != ref_dir)

# learning configuration #
num_train_samples = 100
num_folds = 5
# make sure num_train_samples divides by num folds!
assert (num_train_samples % num_folds == 0)
features = ['t_ref_count', 't_alt_count', 't_lod_fstar', 'tumor_f', 'dbsnp_af', 'dbsnp', 'esp_af',
            'esp', 'thousand_af', 'thousand', 'gnomad_af', 'gnomad',
            'classification_to_remove', 'log_like_RNA', 'log_like_DNA', 'pon_RNA_1', 'pon_RNA_2',
            'pon_RNA_3', 'pon_RNA_4', 'pon_RNA_5', 'pon_RNA_6', 'pon_RNA_7', 'pon_RNA_8',
            'pon_DNA_1', 'pon_DNA_2', 'pon_DNA_3', 'pon_DNA_4', 'pon_DNA_5', 'pon_DNA_6',
            'pon_DNA_7', 'pon_DNA_8']

# environment configuration #
tools = "/Local/md_keren/anaconda3/bin/"
