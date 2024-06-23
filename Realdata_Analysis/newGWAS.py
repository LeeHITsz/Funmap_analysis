import numpy as np
import pandas as pd


# file = 'HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results'
# file = 'LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results'
# file = 'logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results'
# file = 'TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results'
# file = 'without_UKB_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results'
# file = 'without_UKB_LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results'
# file = 'without_UKB_logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results'
file = 'without_UKB_TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results'

data = pd.read_csv("F:/SuSiE/realdata/newGWAS/"+file,
                   usecols=['rsID', 'CHROM', 'pvalue_neg_log10', 'pvalue_neg_log10_GC'], sep='\t')

label = data.loc[data['pvalue_neg_log10'] > -np.log10(5e-8)].values[:, 0]
label_gc = data.loc[data['pvalue_neg_log10_GC'] > -np.log10(5e-8)].values[:, 0]

np.savetxt("newGWAS/Cho_without_UKB_newGWAS.txt", label, fmt='%s')
np.savetxt("newGWAS/Cho_without_UKB_newGWAS_gc.txt", label_gc, fmt='%s')
