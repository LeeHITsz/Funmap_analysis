import pandas as pd


info = pd.read_csv("zscore/variants.tsv.gz", sep='\t')
data = pd.read_csv("zscore/HDL.tsv.gz", sep='\t')
# data = pd.read_csv("zscore/LDL.tsv.gz", sep='\t')
# data = pd.read_csv("zscore/Glucose.tsv.gz", sep='\t')
# data = pd.read_csv("zscore/Cholesterol.tsv.gz", sep='\t')

info = info[['variant', 'chr', 'pos', 'ref', 'alt', 'rsid']]
data = data[['variant', 'minor_allele', 'tstat', 'pval']]

info = info.loc[info['ref'].map(lambda x: len(x) == 1) & info['alt'].map(lambda x: len(x) == 1)]
data = pd.merge(data, info, how='inner', on='variant')

data_pval = data[data["pval"] < 5e-8]

data.to_csv("zscore/HDL.csv.gz", index=False)
data_pval.to_csv("zscore/HDL_pval.csv.gz", index=False)
# data.to_csv("zscore/LDL.csv.gz", index=False)
# data_pval.to_csv("zscore/LDL_pval.csv.gz", index=False)
# data.to_csv("zscore/Glucose.csv.gz", index=False)
# data_pval.to_csv("zscore/Glucose_pval.csv.gz", index=False)
# data.to_csv("zscore/Cholesterol.csv.gz", index=False)
# data_pval.to_csv("zscore/Cholesterol_pval.csv.gz", index=False)
