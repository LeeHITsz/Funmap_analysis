import numpy as np
import pandas as pd
import scipy.sparse as sparse


def load_ld_npz(ld_prefix):
    # load the SNPs metadata
    gz_file = '%s.gz' % (ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid': 'SNP', 'chromosome': 'CHR', 'position': 'BP', 'allele1': 'A1', 'allele2': 'A2'},
                      inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps[
        'A1'] + '.' + df_ld_snps['A2']

    # load the LD matrix
    npz_file = '%s.npz' % (ld_prefix)
    try:
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s' % (npz_file))

    # create df_R and return it
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps


chr_prefix = "ld/chr10_102000001_105000001"
LD, LD_info = load_ld_npz(chr_prefix)

LD_info = LD_info.loc[LD_info['A1'].map(lambda x: len(x) == 1) & LD_info['A2'].map(lambda x: len(x) == 1)]
index = LD_info.index
LD = LD.loc[index][index]
LD_info = LD_info[['SNP', 'A1']]
np.savetxt(chr_prefix +"_ld.txt", LD.values, fmt='%.6f')

