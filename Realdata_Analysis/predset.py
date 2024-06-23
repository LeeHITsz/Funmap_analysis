import os
import numpy as np
import pandas as pd


trait = 'TC'
set_total = pd.read_csv('predset/'+trait+'_predset.csv')

snps_name = np.array([])
for i in range(1, 22+1):

    sets = set_total.loc[set_total['Chromosome'] == i]
    snps = pd.read_csv('zscore/chr'+str(i)+'/zscore/HDL.csv.gz')[['pos', 'rsid']]
    sets = sets.rename(columns={'Position': 'pos'})
    sets = pd.merge(sets, snps, how='inner', on='pos')
    snps_name = np.append(snps_name, sets['rsid'].values)

np.savetxt("predset/"+trait+"_predset.txt", snps_name, fmt='%s')

# def get_all_files(directory):
#
#     file_paths = []
#
#     for root, dirs, files in os.walk(directory):
#         for file in files:
#             file_path = os.path.join(root, file)
#             file_paths.append(file_path)
#
#     return file_paths
#
#
# directory_path = "predset/logTG"
# load_list = get_all_files(directory_path)
#
# sets = pd.DataFrame(columns=['Chromosome', 'Position'])
# for load in load_list:
#     sets_tmp = pd.read_csv(load, sep='\t').iloc[:, 0:2]
#     sets = pd.concat((sets, sets_tmp))
# sets.to_csv('predset/TG_predset.csv', index=False)
