import tarfile

genotype_compressed = '../data/LDREF.tar.bz2'
phenotype_data = '../data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz'

def extract_data():
    print('Extracting data...')
    tar = tarfile.open(genotype_compressed, mode='r:bz2')
    tar.extractall('../data')
    tar.close()
    print('Done')