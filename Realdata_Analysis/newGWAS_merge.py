import numpy as np
import pandas as pd


threshold = 0.95
trait_list = ['HDL', 'LDL', 'TG', 'TC']
method_list = ['SuSiE', 'CARMA', 'PAINTOR', 'Funmap', 'CARMA_anno', 'PAINTOR_anno']


for trait in trait_list:

    # _gc, _without_UKB, _without_UKB_gc,
    GWAS = np.loadtxt(('newGWAS/'+trait+'_MR.txt'), dtype=str)
    GWAS = pd.Series(GWAS, name='rsID')

    for method in ['Funmap', 'CARMA', 'PAINTOR']:

        pip = pd.read_csv('newGWAS/'+method+'_'+trait.replace('TC', 'Cho')+'.csv', header=None)
        pip.columns = ['rsID', 'PIP']
        pip = pip.loc[pip['PIP'] > threshold]
        merge = pd.merge(pip, GWAS, how='inner', on='rsID')

        print(method, ',', trait, ':', merge.shape[0], '/', pip.shape[0]-merge.shape[0], '=', merge.shape[0]/pip.shape[0])


# trait = 'TG'
# snp_susie = pd.read_csv('newGWAS/' + 'SuSiE_' + trait + '80.csv', header=None)
# snp_funmap = pd.read_csv('newGWAS/' + 'Funmap_' + trait + '80.csv', header=None)
# snp_funmap[~np.isin(snp_funmap.iloc[:, 0], snp_susie.iloc[:, 0])].to_csv('newGWAS/Funmap_' + trait + '.csv',
#                                                                          index=False, header=None)
#
# snp_carma = pd.read_csv('newGWAS/' + 'CARMA_' + trait + '80.csv', header=None)
# snp_carma_anno = pd.read_csv('newGWAS/' + 'CARMA_anno_' + trait + '80.csv', header=None)
# snp_carma_anno[~np.isin(snp_carma_anno.iloc[:, 0], snp_carma.iloc[:, 0])].to_csv('newGWAS/CARMA_' + trait + '.csv',
#                                                                                  index=False, header=None)
#
# snp_paintor = pd.read_csv('newGWAS/' + 'PAINTOR_' + trait + '80.csv', header=None)
# snp_paintor_anno = pd.read_csv('newGWAS/' + 'PAINTOR_anno_' + trait + '80.csv', header=None)
# snp_paintor_anno[~np.isin(snp_paintor_anno.iloc[:, 0], snp_paintor.iloc[:, 0])].to_csv(
#     'newGWAS/PAINTOR_' + trait + '.csv', index=False, header=None)
