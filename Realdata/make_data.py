import os
import time
import subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sparse
from functools import reduce


Chr = '22'

anno = pd.read_csv('chr'+Chr+'/annotation/baselineLF.'+Chr+'.annot', sep='\t')
anno.drop(anno[['CHR', 'BP', 'CM']], axis=1, inplace=True)
anno.drop_duplicates(subset=['SNP'], keep='first', inplace=True)

ld_list = np.loadtxt('chr'+Chr+'/LD/LD_list.txt', dtype=int)
ld_list = ld_list[:]

z_HDL = pd.read_csv('chr'+Chr+'/zscore/HDL.csv.gz').drop_duplicates(subset=['rsid'], keep='first')
z_LDL = pd.read_csv('chr'+Chr+'/zscore/LDL.csv.gz').drop_duplicates(subset=['rsid'], keep='first')
z_TG = pd.read_csv('chr'+Chr+'/zscore/TG.csv.gz').drop_duplicates(subset=['rsid'], keep='first')
z_Cho = pd.read_csv('chr'+Chr+'/zscore/Cholesterol.csv.gz').drop_duplicates(subset=['rsid'], keep='first')

z_HDL_pval = pd.read_csv('chr'+Chr+'/zscore/HDL_pval.csv.gz').drop_duplicates(subset=['rsid'], keep='first')
z_LDL_pval = pd.read_csv('chr'+Chr+'/zscore/LDL_pval.csv.gz').drop_duplicates(subset=['rsid'], keep='first')
z_TG_pval = pd.read_csv('chr'+Chr+'/zscore/TG_pval.csv.gz').drop_duplicates(subset=['rsid'], keep='first')
z_Cho_pval = pd.read_csv('chr'+Chr+'/zscore/Cholesterol_pval.csv.gz').drop_duplicates(subset=['rsid'], keep='first')


def get_ldname(chr, pos):

    return 'chr'+str(chr)+'_'+str(int(pos*1e6+1))+'_'+str(int((pos+3)*1e6+1))


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


for Pos in ld_list:

    if not os.path.exists('chr' + Chr + '/LD/' + get_ldname(Chr, Pos) + '.gz'):

        cmd = 'aws s3api get-object --bucket broad-alkesgroup-ukbb-ld ' + \
              '--key UKBB_LD/' + get_ldname(Chr, Pos) + '.gz ' \
              'chr' + Chr + '/LD/' + get_ldname(Chr, Pos) + '.gz'

        while True:
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, err = process.communicate()
            if process.returncode == 0:
                break
            else:
                print("Download failed, retrying...")
                time.sleep(5)

    if not os.path.exists('chr' + Chr + '/LD/' + get_ldname(Chr, Pos) + '.npz'):

        cmd = 'aws s3api get-object --bucket broad-alkesgroup-ukbb-ld ' + \
              '--key UKBB_LD/' + get_ldname(Chr, Pos) + '.npz ' \
              'chr' + Chr + '/LD/' + get_ldname(Chr, Pos) + '.npz'

        while True:
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, err = process.communicate()
            if process.returncode == 0:
                break
            else:
                print("Download failed, retrying...")
                time.sleep(5)

    LD, LD_info = load_ld_npz('chr' + Chr + '/LD/' + get_ldname(Chr, Pos))
    LD_info.drop_duplicates(subset=['SNP'], keep='first', inplace=True)
    LD_info = LD_info.loc[(LD_info['BP'] < ((Pos+2)*1e6)) & (LD_info['BP'] > ((Pos+1)*1e6))]
    LD_info = LD_info.loc[LD_info['A1'].map(lambda x: len(x) == 1) & LD_info['A2'].map(lambda x: len(x) == 1)]

    SNP_list = pd.Series()

    for z, z_pval, trait in zip([z_HDL, z_LDL, z_TG, z_Cho], [z_HDL_pval, z_LDL_pval, z_TG_pval, z_Cho_pval],
                                ['HDL', 'LDL', 'TG', 'Cho']):

        if (Pos == 0) and not any([1 <= x <= 2e6 + 1 for x in z_pval['pos'].values]):

            continue

        elif not any([(Pos + 1) * 1e6 + 1 <= x <= (Pos + 2) * 1e6 + 1 for x in z_pval['pos'].values]):

            continue

        else:

            if SNP_list.empty:

                SNP_list = pd.Series(reduce(np.intersect1d, (LD_info['SNP'], anno['SNP'], z['rsid'])), name='SNP')
                LD_info = LD_info[['SNP', 'A1']]
                SNP_index = LD_info[LD_info['SNP'].isin(SNP_list.values)].index
                LD_info = pd.merge(LD_info, SNP_list, how='inner', on='SNP')

                LD = LD.loc[SNP_index][SNP_index]
                np.savetxt('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '.ld', LD.values, fmt='%.3f')

                anno_pos = pd.merge(anno, SNP_list, how='inner', on='SNP')
                anno_pos.drop(anno[['SNP']], axis=1, inplace=True)
                anno_pos.to_csv('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '.annotations', index=False, sep=' ')

            # z.rename(columns={'rsid': 'SNP'}, inplace=True)
            z = pd.merge(z, LD_info, how='inner', left_on='rsid', right_on='SNP')
            inv_index = (z['A1'] != z['alt'])
            z.loc[inv_index, 'tstat'] = -z.loc[inv_index, 'tstat']

            z_score = pd.DataFrame({'SNP': z['rsid'].values, 'z': z['tstat'].values})
            z_score.to_csv('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '_' + trait + '.txt',
                           header=None, index=False, sep='\t', float_format='%.6f')


