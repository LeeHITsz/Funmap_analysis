import pandas as pd

# data = pd.read_csv("zscore/HDL.csv.gz")
# data = pd.read_csv("zscore/LDL.csv.gz")
# data = pd.read_csv("zscore/Glucose.csv.gz")
# data = pd.read_csv("zscore/Cholesterol.csv.gz")
data = pd.read_csv("zscore/TG.csv.gz")

# data_pval = pd.read_csv("zscore/HDL_pval.csv.gz")
# data_pval = pd.read_csv("zscore/LDL_pval.csv.gz")
# data_pval = pd.read_csv("zscore/Glucose_pval.csv.gz")
# data_pval = pd.read_csv("zscore/Cholesterol_pval.csv.gz")
data_pval = pd.read_csv("zscore/TG_pval.csv.gz")

for i in range(1, 22+1):

    data_tmp = data.loc[data['chr'].map(lambda x: str(x) == str(i))]
    data_pval_tmp = data_pval.loc[data_pval['chr'].map(lambda x: str(x) == str(i))]
    data_tmp.to_csv("zscore/chr"+str(i)+"/zscore/TG.csv.gz", index=False)
    data_pval_tmp.to_csv("zscore/chr"+str(i)+"/zscore/TG_pval.csv.gz", index=False)

    if not data_pval_tmp['pos'].empty:
        print(i, ": ", max(data_pval_tmp['pos'].values))
    else:
        print(i, ": ", 0)
