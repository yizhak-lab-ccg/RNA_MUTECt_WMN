---------------
RNA-MuTect-WMN
---------------

*****************************************
The pipeline uses Python v3.9.1
*****************************************

No installation is required.

INPUT:
The pipeline includes the following input files that should be placed in a folder named 'input_folder' (see example input files in 'input_folder'):
1. 5 text files for train sets 1-5, containing variants in the rows and all features in the column. Last column is the label (somatic=1; no somatic=0). 
2. A similar text file containing all variants in the test set.

OUTPUT:
1. 5 trained models
2. Precision and recall values for train and test groups
* Running time is dependent on the number of variants in the dataset



