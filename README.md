# Genetic Risk Prediction with cis-eQTLs and Polygenic Risk Scores
## Garvey Li

*README and instructions still WIP*


**Instructions**

1. The majority of the running code is located in `src/single_gene_prs.ipynb`. Install all the necessary requirements to run this notebook using `pip install -r requirements.txt` in the root directory of the repository. This can also be done in `src/single_gene_prs.ipynb` by uncommenting the line of code in the first cell and running the cell. 

2. Run the 2nd cell to import the required packages as well as load the data required for my tools to run.

3. Run the 3rd cell to ensure that the `plink2` software has been installed and is running correctly. More specific instructions on what to do are written in the notebook.

4. Load the datasets we are going to look at. There are a total of 4 dataframes
    * alleles: A dataframe containing SNP metadata, such as chromosome, SNP id, position, and the reference and effect alleles. 
    * samples: A dataframe containing metadata of the individuals whose data was used in the 1000 Genomes Project.
    * genotypes: A dataframe containing the genetic variation measured for each individual across 1,190,321 SNPs
    * P (phenotypes): A dataframe containing the levels of gene expression for each individual across 23,722 genes. 

    

5. To run an eQTL analysis on any gene in cell 5, simply pick a gene from the `P` dataframe and store its `Chr` and `TargetID` values in the variables `chr_to_explore` and `gene_to_explore`, respectively.

6. To run an eQTL across all genes in a chromosome, uncomment the 2 lines of code in cell 6, and run the cell.

7. To create a histogram of the polygenic risk scores of the gene selected in step 4, run cell 7.