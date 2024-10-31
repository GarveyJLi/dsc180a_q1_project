import pandas_plink as pp
import statsmodels.api as sm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from functools import reduce
import os
import tarfile
import gzip

print('Loading data...')
phenotype_data = '../data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz'
# alleles = .bim (snp. came from vcf file (txt file))
# samples = .fam (family id)
# genotypes = .bed (genotypes). Rows and columns defined by alleles and samples
(alleles, samples, genotypes) = pp.read_plink("../data/LDREF/1000G.EUR.*",
                             verbose=False)
genotypes = pd.DataFrame(genotypes.compute())
P = pd.read_csv(phenotype_data, sep='\t', compression='gzip')
print('Complete')

def get_chr(gene):
    target_gene_info = P[P['TargetID'] == gene].iloc[0]
    target_chr = str(target_gene_info['Chr'])
    return target_chr

def generate_pos_file(gene, out_dir='', window=500000):
    """
    Generates a position file from a gene to be used in get_snps() (or plink2 -extract)
    """
    target_gene_info = P[P['TargetID'] == gene].iloc[0]
    target_chr = str(target_gene_info['Chr'])
    target_start = target_gene_info['Coord'] - window
    target_end = target_gene_info['Coord'] + window
    coord_str = f'{target_chr}\t{target_start}\t{target_end}\t{gene}\t0\t+'
    file_name = f'{out_dir}{gene}_coord.txt'
    with open(file_name, "w") as text_file:
        text_file.write(coord_str)
    return file_name

def get_snps(gene, pos_file='', out_dir='', delete_pos=True):
    chr = get_chr(gene)
    if pos_file == '':
        pos_file = generate_pos_file(gene, out_dir=out_dir)
    plink_cmd = f'./plink2 --bfile ../data/LDREF/1000G.EUR.{chr} --extract bed1 {pos_file} --out {out_dir}1000G.EUR.{chr}.{gene} --make-bed'
    os.system(plink_cmd)
    if delete_pos:
        os.remove(pos_file)
    return plink_cmd

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
    target_chr_data = phenotypes[phenotypes['Chr']==chr]
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

    #return pd.concat(all_analyses).merge(pd.merge(alleles[['snp', 'chrom', 'pos', 'a0', 'a1', 'i']], effect_allele_freq, how='left', on='snp'), how='left', on='snp')

