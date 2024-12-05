# Genetic Risk Prediction with cis-eQTLs and Polygenic Risk Scores
## Garvey Li


## Welcome!
Hi there! Welcome to my genetic risk prediction project. This repository contains the code I used to analyze genetic associations of cis-eQTLs across the human genome, as well as creating genetic risk models to predict gene expression levels and traits.

This repository contains three notebooks corresponding to the 3 analyses done in my paper, which is linked [here](genetic_risk_prediction_with_prs.pdf). The sections are as follows:

1. Genome-wide cis-eQTL analysis
2. Single gene PRS model training and evaluation using cis-eQTL analyses on the 1000 Genomes dataset
3. PRS model training using GWAS summary statistics pertaining to individuals' heights.


## Instructions

1. The data extraction and package installations required for the code to run can be done in any of the 3 notebooks located in the `src\` folder by uncommenting 
`!pip install -r ../requirements.txt` in the first python cell and `extract_data()` in the second python cell, then running both cells. This only needs to be done once, so no need to do it in more than 1 notebook. 

2. After the packages have been installed and the data has been extracted, there is another step to take to ensure the tools written for the analyses work properly. The code requires a piece of software name `plink` and `plink2`. By default, it may not have permissions that allow it to work, so the following lines of code in cell 3 must be uncommented and run:

    `!chmod 700 plink2`

    `!chmod 700 plink`

3. Once you've made sure the required packages have been installed, the data has been extracted, and plink is able to run, you should be able to run the code in the notebooks.

### Part 1: Genome-wide cis-eQTL analysis

This section of my analysis contains results from the notebook located at `src\genome_wide_ciseqtl_analysis.ipynb`. 

The code in cell 5 originally parsed through each of the 23722 genes available in the dataset across 22 chromosomes in order to conduct cis-eQTL analyses on each of them. However, the code for this has been commented out and the results of this have been saved to `src\genome_wide_ciseqtl_analysis\`, as it takes quite a while to run. 

The results can be loaded in the 6th python cell, titled `Results from the genome wide cis-eQTL analysis`, and may take a few seconds before completion, as the result files are so large. If this code fails to run at all, it is likely that your machine has run out of memory, so be wary.

### Part 2: Single gene PRS model training and evaluation using cis-eQTL analyses on the 1000 Genomes dataset

This section of my analysis contains results from the notebook located at `src\single_gene_prs.ipynb`. 

The code in this notebook generates PRS models for various genes related to cardiovascular diseases, type 1 diabetes, and type 2 diabetes (in that order). They are evaluated against their respective null models and compared to the corresponding gene expression distributions. 

The code required to generate 10 different PRS distributions requires quite a bit of memory, so once again, be wary.

### Part 3: PRS model training using GWAS summary statistics pertaining to individuals' heights.

This section of my analysis contains results from the notebook located at `src\height_prs.ipynb`. 

This section explores the use of genome-wide association study (GWAS) summary statistics to train PRS models predicting individuals' heights. Although the summary statistics do not contain any data that allow for evaluation (e.g. gene expression levels), they can be loosely evaluated using individual's samples with known traits. A `.vcf` file created from 23andMe data can be uploaded to the `data\` folder and compared against the PRS model for height. To do so, replace the vcf filepath `vcf_fp = '../data/f_eur_sample.vcf'` in cell 4 with the name of your `.vcf` file, and you can see where you lie on the PRS distribution. 

