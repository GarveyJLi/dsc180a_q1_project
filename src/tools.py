import pandas_plink as pp
import statsmodels.api as sm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from functools import reduce
import os
import tarfile
import gzip


genotype_compressed = '../data/LDREF.tar.bz2'
phenotype_data = '../data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz'

def extract_data():
    print('Extracting data...')
    tar = tarfile.open(genotype_compressed, mode='r:bz2')
    tar.extractall('../data')
    tar.close()
    print('Done')



# alleles = .bim (snp. came from vcf file (txt file))
# samples = .fam (family id)
# genotypes = .bed (genotypes). Rows and columns defined by alleles and samples
(alleles, samples, genotypes) = pp.read_plink("../data/LDREF/1000G.EUR.*",
                             verbose=False)
genotypes = pd.DataFrame(genotypes.compute())

#allele_freq_df = genotypes.copy(deep=True)
#allele_freq_df['snp'] = alleles['snp']
#allele_freq_df['freq_0'] = genotypes.apply(lambda x: x == 0).sum(axis=1) / 489
#allele_freq_df['eaf'] = genotypes.apply(lambda x: x != 0).sum(axis=1) / 489
#effect_allele_freq = allele_freq_df[['snp', 'eaf']]
P = pd.read_csv(phenotype_data, sep='\t', compression='gzip')
P

def cis_eQTL_analysis(chr, gene, alleles, samples, genotypes, phenotypes, window=500000, mult_gene=False):
    target_gene_data = phenotypes[phenotypes['Gene_Symbol']==gene].reset_index(drop=True)
    target_gene_data['start'] = target_gene_data['Coord'] - window
    target_gene_data['end'] = target_gene_data['Coord'] + window
    target_gene_alleles = alleles[(alleles['chrom'] == chr) & (alleles['pos'] <= target_gene_data['end'].values[0]) & (alleles['pos'] >= target_gene_data['start'].values[0])]
    snp_ids = dict(zip(target_gene_alleles.index, target_gene_alleles['snp'].values))
    target_gene_genotypes = genotypes.iloc[target_gene_alleles.index].T.rename(columns=snp_ids)
    sample_ids = [t for t in target_gene_data.columns if t not in ['TargetID', 'Gene_Symbol', 'Chr', 'start', 'end']]
    full_target_gene_data = target_gene_genotypes.merge(samples, how='left', left_index=True, right_index=True).merge(target_gene_data[sample_ids].T.rename(columns=snp_ids), how='left', left_on='fid', right_index=True).rename(columns={0: 'expression'}).drop(columns=['iid', 'father', 'mother', 'gender', 'trait', 'i'])
    full_target_gene_data = full_target_gene_data[~full_target_gene_data['expression'].isna()]

    y = full_target_gene_data['expression'].to_numpy()
    all_snps = full_target_gene_data.columns[:-2]
    snp_regs = {}
    snp_data = {}
    for snp in all_snps:
        temp_x = full_target_gene_data[snp].values 
        temp_x = sm.add_constant(temp_x, has_constant='add')
       
        model = sm.OLS(y,temp_x)
        results = model.fit()
        snp_regs[snp] = results
        snp_data[snp] = {
            'p_val': results.pvalues[1],    
            'beta': results.params[1],
            'se': results.tvalues[1],
        } 

    target_gene_linreg = pd.DataFrame.from_dict(snp_data, orient='index').reset_index().rename(columns={'index': 'snp'})
    target_gene_linreg['gene'] = gene
    target_gene_linreg['gene_pos'] = target_gene_data['Coord'].values[0]
    target_gene_linreg['n'] = samples.shape[0]
    
    if mult_gene:
        return target_gene_linreg
    return target_gene_linreg.merge(alleles[['snp', 'chrom', 'pos', 'a0', 'a1', 'i']], how='left', on='snp')




def chr_analysis(chr, alleles, samples, genotypes, phenotypes, window=500000):
    target_chr_data = P[P['Chr']==chr]
    all_genes = target_chr_data['Gene_Symbol'].unique()
    all_analyses=[]
    print(f'Starting cis-eQTL analyses for chromosome {chr}')

    gene_counter = 0
    total_genes = len(all_genes)
    for g in all_genes:
        if gene_counter%100==0:
            print(f'{gene_counter}/{total_genes} analyses completed')
        temp_df = cis_eQTL_analysis(chr, g, alleles, samples, genotypes, phenotypes, window, True)
        all_analyses.append(temp_df.copy())
        gene_counter += 1
    print(f'Compiling results')
    return pd.concat(all_analyses).merge(alleles[['snp', 'chrom', 'pos', 'a0', 'a1', 'i']], how='left', on='snp')

