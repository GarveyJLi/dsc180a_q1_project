import pandas_plink as pp
import statsmodels.api as sm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from functools import reduce
import random
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

def r_squared(df, x1, x2):
    temp_r = df[x1].corr(df[x2])
    if np.isnan(temp_r):
        temp_r = 0
    return temp_r ** 2


def get_chrom(gene):
    target_gene_info = P[P['TargetID'] == gene].iloc[0]
    target_chrom = str(target_gene_info['Chr'])
    return target_chrom

def generate_pos_file(gene, out_dir='', window=500000):
    """
    Generates a position file from a gene to be used in get_snps() (or plink2 -extract)
    """
    target_gene_info = P[P['TargetID'] == gene].iloc[0]
    target_chrom = str(target_gene_info['Chr'])
    target_start = target_gene_info['Coord'] - window
    target_end = target_gene_info['Coord'] + window
    coord_str = f'{target_chrom}\t{target_start}\t{target_end}\t{gene}\t0\t+'
    file_name = f'{out_dir}{gene}_coord.txt'
    with open(file_name, "w") as text_file:
        text_file.write(coord_str)
    return file_name

def get_snps(gene, pos_file='', out_dir='', delete_pos=True):
    chrom = get_chrom(gene)
    if pos_file == '':
        pos_file = generate_pos_file(gene, out_dir=out_dir)
    plink_cmd = f'{os.getcwd()}/plink2 --bfile ../data/LDREF/1000G.EUR.{chrom} --extract bed1 {pos_file} --out {out_dir}1000G.EUR.{chrom}.{gene} --make-bed'
    os.system(plink_cmd)
    if delete_pos:
        os.remove(pos_file)
    return plink_cmd

def cis_eQTL_analysis(chrom, gene, alleles, samples, genotypes, phenotypes, window=500000, mult_gene=False):
    target_gene_data = phenotypes[phenotypes['Gene_Symbol']==gene].reset_index(drop=True)
    target_gene_data['start'] = target_gene_data['Coord'] - window
    target_gene_data['end'] = target_gene_data['Coord'] + window
    target_gene_alleles = alleles[(alleles['chrom'] == chrom) & (alleles['pos'] <= target_gene_data['end'].values[0]) & (alleles['pos'] >= target_gene_data['start'].values[0])]
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
            'p': results.pvalues[1],    
            'beta': results.params[1],
            'se': results.tvalues[1],
        } 

    target_gene_linreg = pd.DataFrame.from_dict(snp_data, orient='index').reset_index().rename(columns={'index': 'snp'})
    target_gene_linreg['gene'] = gene
    target_gene_linreg['gene_pos'] = target_gene_data['Coord'].values[0]
    target_gene_linreg['n'] = samples.shape[0]
    
    if mult_gene:
        return target_gene_linreg
    return target_gene_linreg.merge(alleles[['snp', 'chrom', 'pos', 'a0', 'a1', 'i']], how='left', on='snp').rename(columns={'chrom': 'chr', 'pos': 'bp'})[['snp', 'chr', 'bp', 'p', 'se', 'beta','a0', 'a1', 'i']]


def chr_analysis(chrom, alleles, samples, genotypes, phenotypes, window=500000):
    target_chrom_data = phenotypes[phenotypes['Chr']==chrom]
    all_genes = target_chrom_data['Gene_Symbol'].unique()
    all_analyses=[]
    print(f'Starting cis-eQTL analyses for chromosome {chr}')

    gene_counter = 0
    total_genes = len(all_genes)
    for g in all_genes:
        if gene_counter%100==0:
            print(f'{gene_counter}/{total_genes} analyses completed')
        temp_df = cis_eQTL_analysis(chrom, g, alleles, samples, genotypes, phenotypes, window, True)
        all_analyses.append(temp_df.copy())
        gene_counter += 1
    print(f'Compiling results')
    return pd.concat(all_analyses).merge(alleles[['snp', 'chrom', 'pos', 'a0', 'a1', 'i']], how='left', on='snp').rename(columns={'chrom': 'chr', 'pos': 'bp'})[['snp', 'chr', 'bp', 'p', 'se', 'beta', 'a0', 'a1', 'i']]

def clump(gene, eqtl_analysis_fp, p_val, r2, out_dir=''):
    # Get all snps within 500 kB of the target gene
    get_snps(gene, out_dir=f'{os.getcwd()}/temp_files/', delete_pos=True)
    # Clump snps
    to_clump = f'{os.getcwd()}/temp_files/1000G.EUR.{get_chrom(gene)}.{gene}'
    clump_cmd = f'{os.getcwd()}/plink2 --bfile {to_clump} --clump-p1 {p_val} --clump-r2 {r2} --clump-kb 250 --clump {eqtl_analysis_fp} --clump-snp-field snp --clump-field p --out {to_clump}'
    os.system(clump_cmd)
    #extract SNPs from clump output
    out_snps = f'{os.getcwd()}/temp_files/PRS.SNPs.{gene}'
    extract_clump_cmd = "awk 'NR!=1{print $3}' " + f'{to_clump}.clumps > {out_snps}'
    os.system(extract_clump_cmd)
    return pd.read_csv(out_snps, header=None).rename(columns={0:'snp'})


def filter_snps(all_snps, eqtl_analysis, alleles, genotypes):
    filtered_eqtl_analysis = eqtl_analysis[eqtl_analysis['snp'].isin(all_snps['snp'])]
    filtered_alleles = alleles[alleles['snp'].isin(all_snps['snp'])]
    filtered_genotypes = genotypes.loc[filtered_alleles['i']]
    return filtered_eqtl_analysis, filtered_alleles, filtered_genotypes

def prs_gene_pipeline(analysis_df, samples, genotypes):
    #analysis_df = analysis_df[(analysis_df['p'] <= p_val)]
    prs = {}
    for iid in samples['i']:
        numerator = 0
        denominator = 2 * analysis_df.shape[0]
        for snp in analysis_df['i']:
            temp_snp = analysis_df[analysis_df['i']==snp]
            temp_beta = temp_snp['beta'].values[0]
            temp_ea_count = genotypes[iid].loc[snp]
            numerator += temp_beta * temp_ea_count
        temp_fid = samples[samples['i']==iid]['fid'].values[0]
        
        prs[temp_fid] = numerator / (denominator + 1e-10)
    return prs
    
    
def generate_prs(gene, alleles, samples, genotypes, phenotypes, p_val, r2, random_seed=0):
    np.random.seed(random_seed)
    # Sampling population
    known_expressions = phenotypes[phenotypes['TargetID']==gene]
    known_fids = list(known_expressions.columns[4:].values)
    random.shuffle(known_fids)
    
    # Train valid test splits. 80 10 10
    num_analysis = int(0.8 * len(known_fids))
    num_valid = (len(known_fids) - num_analysis) // 2
    num_test = len(known_fids) - num_valid - num_analysis
    analysis_fids = known_fids[:num_analysis]
    valid_fids = known_fids[num_analysis:num_analysis+num_valid]
    test_fids = known_fids[num_analysis+num_valid:]

    analysis_samples = samples[samples['fid'].isin(analysis_fids)]
    analysis_genotypes = genotypes[analysis_samples['i'].values]
    analysis_phenotypes = phenotypes[phenotypes['TargetID']==gene][['TargetID', 'Gene_Symbol', 'Chr', 'Coord'] + list(analysis_fids)]

    valid_samples = samples[samples['fid'].isin(valid_fids)]
    valid_genotypes = genotypes[valid_samples['i'].values]
    valid_phenotypes = phenotypes[phenotypes['TargetID']==gene][['TargetID', 'Gene_Symbol', 'Chr', 'Coord'] + list(valid_fids)]


    test_samples = samples[samples['fid'].isin(test_fids)]
    test_genotypes = genotypes[test_samples['i'].values]
    test_phenotypes = phenotypes[phenotypes['TargetID']==gene][['TargetID', 'Gene_Symbol', 'Chr', 'Coord'] + list(test_fids)]
    
    # Generate initial eQTL analysis results
    target_gene_analysis = cis_eQTL_analysis(get_chrom(gene), gene, \
                            alleles, analysis_samples, analysis_genotypes, \
                                analysis_phenotypes)
    # Write results to csv for clumping
    analysis_out = os.getcwd() + f'/temp_files/{gene+"_eqtl_analysis.txt"}'
    target_gene_analysis.to_csv(analysis_out, index=False, sep='\t')
    
    # Clumping
    clumped_snps_lst = clump(gene, analysis_out, p_val, r2, out_dir=f'{os.getcwd()}/temp_files')
    clumped_eqtl, clumped_alleles, clumped_genotypes = filter_snps(clumped_snps_lst, target_gene_analysis, alleles, genotypes)

    prs_scores = prs_gene_pipeline(clumped_eqtl, samples, clumped_genotypes)
    return prs_scores                            



    
