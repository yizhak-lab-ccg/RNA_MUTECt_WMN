# RNA-Mutect-WMN
This pipeline implements the method described in 
[Estimating tumor mutational burden from RNA-sequencing without a matched-normal sample](https://www.nature.com/articles/s41467-022-30753-2), 
and **should be used after running [RNA_MuTect](https://github.com/broadinstitute/RNA_MUTECT_1.0-1).**
This pipeline runs on a Linux machine only.

## Requirements
1. python3 packages:
   * pandas (1.1.5+)
   * NumPy (1.19.4+)
   * scikit-learn (0.23.2+)
   * matplotlib (3.3.3+)
2. [CAPY](https://github.com/getzlab/CApy/tree/master/capy) python package (0.1+)
3. [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) as part of the [gatk](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) package (4.2.6.1+)
4. Samtools: 
   * [bgzip](http://www.htslib.org/doc/bgzip.html) (1.11+)
   * [bcftools](https://samtools.github.io/bcftools/bcftools.html) (1.8+)
   * [tabix](http://www.htslib.org/doc/tabix.html) (1.11+)
5. ~300 GB space: the 'resource' folder will be around 230 GB, and more space will be required (depending on the number of samples).

## Input files and directory tree
Directory names can be changed in the [configuration](#configuration) file
```
Data/
    'cancer_dir'/                               #project-specific
         input/
            call_stats/
            maf/
    resource
        BCF_tools_dbs/
            merged.vcf.gz                        #ESP db
        pon/
            'RNA_binary'                              
            'DNA_binray'
        reference/
            'reference.fasta'
            'reference.fasta.fai'
            'reference.dict'          
```

## Configuration
The `config.py` file should be configured by the user.
1. [**Directory configuration**](#input-files-and-directory-tree):
   1. 'cancer_type' is the name of the project-specific directory.
   2. other directories and file names can be changed using this file if desired.
2. **Learning configuration**: in this section, you can play with the learning parameters and features.
3. **Environment configuration** is used to configure some tools' locations.
    1. `tools` is the location of the [samtools and GATK](#requirements) binaries

## Running instructions

### Inputs preparation
1. As mentioned before, the input of this tool is the output of [RNA-MuTect](https://github.com/broadinstitute/RNA_MUTECT_1.0-1). 
A cloud implementation can be found in [Terra](https://app.terra.bio/#workspaces/broad-firecloud-gtex/RNA_MuTect).
    * In order to run RNA-MuTect any normal sample can be used and it does not require the matched-normal sample.
2. Details for location of PoN files are in the [manuscript](https://www.nature.com/articles/s41467-022-30753-2) under 'Data Availability'.
3. The human reference genome hg19 reference files should be used.
4. After downloading the repo, [directory configuration](#Input-files-and-directory-tree) should be done, using the `config.py` file:
   * Under the 'Data' folder:
     * create a 'cancer_dir' folder and configure its name in `config.py`.
   * Under the 'cancer_dir' folder:
     * Create an 'input' folder, and under it a 'maf' and 'call_stats' folders. 
     * Download 'call_stats_capture_paper_v1_3' files (RNA-MuTect output) into 'call_stats' folder.
     * Download 'maf_file_rna_final_paper_v1_3' files (RNA-MuTect output) into 'maf folder'.
   * Under the 'resource' folder:
     * download the pon binary files (DNA & RNA) into the 'pon' folder
     * download the reference files (including .fasta.fai and .dict files) into the 'reference' folder
     * configure downloaded file names in `config.py`.


### Run pipeline
Run `RNA-Mutect-WMN.py`

## Results
When the tool is finished successfully, a 'results' directory will be created under the specified 'cancer_dir'. inside 'results' directory:
1. Train results
   1. mean recall and precision scores
   2. mean recall and precision scores **per sample** + boxplot
2. 'somatics.maf': MAF file of all the variants classified as somatic by the tool. This should be further filtered using RNA-MuTect filtering steps as described in the [paper](https://www.nature.com/articles/s41467-022-30753-2)
