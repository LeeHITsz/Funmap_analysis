import numpy as np
import pandas as pd


pos_max = np.loadtxt("chr_pos_TG.txt")

for i in range(0, 22):

    ld_list = []

    for j in range(int(pos_max[i]//1e6)):

        # for z in ['HDL', 'LDL', 'Glucose', 'Cholesterol']:
        #
        #     pos_list = pd.read_csv('zscore/chr'+str(i+1)+'/zscore/'+z+'_pval.csv.gz')['pos'].values
        #
        #     if j == 0 and any([1 <= x <= 2e6+1 for x in pos_list]):
        #         ld_list.append(j)
        #         break
        #
        #     elif any([(j+1)*1e6+1 <= x <= (j+2)*1e6+1 for x in pos_list]):
        #         ld_list.append(j)
        #         break

        z = 'TG'
        pos_list = pd.read_csv('zscore/chr'+str(i+1)+'/zscore/'+z+'_pval.csv.gz')['pos'].values

        if j == 0 and any([1 <= x <= 2e6+1 for x in pos_list]):
            ld_list.append(j)

        elif any([(j+1)*1e6+1 <= x <= (j+2)*1e6+1 for x in pos_list]):
            ld_list.append(j)

    print("chr"+str(i+1)+": "+str(len(ld_list)))
    np.savetxt('zscore/chr'+str(i+1)+'/LD/LD_list_TG.txt', ld_list, fmt='%d')
