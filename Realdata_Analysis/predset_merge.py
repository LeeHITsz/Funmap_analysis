import numpy as np
import pandas as pd


threshold = 0.95
trait_list = ['HDL', 'LDL', 'TG', 'TC']
method_list = ['SuSiE', 'CARMA', 'PAINTOR', 'Funmap', 'CARMA_anno', 'PAINTOR_anno']


for trait in trait_list:

    # _gc, _without_UKB, _without_UKB_gc,
    sets = np.loadtxt('predset/'+trait+'_predset.txt', dtype=str)
    sets = pd.Series(sets, name='rsID')

    for method in ['Funmap', 'CARMA', 'PAINTOR']:

        pip = pd.read_csv('predset/'+method+'_'+trait+'.csv', header=None)
        pip.columns = ['rsID', 'PIP']
        pip = pip.loc[pip['PIP'] > threshold]
        merge = pd.merge(pip, sets, how='inner', on='rsID')

        print(method, ',', trait, ':', merge.shape[0], '/', pip.shape[0]-merge.shape[0], '=', merge.shape[0]/pip.shape[0])
