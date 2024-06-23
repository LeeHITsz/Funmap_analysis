import os
import time
import numpy as np
import pandas as pd

index = 'Timing'
gene = 'TBX3'  # 'GENE8', 'GENE5', 'TBX3', 'C6orf211', 'MYC'
label = '9'
m = 200

LD = pd.read_csv(index + "/" + gene + '/data/ld.txt', sep='\s+', header=None)
np.savetxt(index + "/" + gene + "/data/LD.txt", LD, fmt='%.3f')

anno = ["X" + str(j) for j in range(1, m+1)]

# cmd1 = "PAINTOR_V3.0/PAINTOR " \
#        "-input " + index + "/" + gene + "/data/input.files " \
#        "-in " + index + "/" + gene + "/data/ " \
#        "-out " + index + "/" + gene + "/PAINTOR_result_no/ " \
#        "-Zhead Zscore " \
#        "-LDname ld " \
#        "-mcmc " \
#        "-annotations Coding"

cmd2 = "PAINTOR_V3.0/PAINTOR " \
       "-input " + index + "/" + gene + "/data/input.files " \
       "-in " + index + "/" + gene + "/data/ " \
       "-out " + index + "/" + gene + "/PAINTOR_result_anno/ " \
       "-Zhead Zscore " \
       "-LDname ld " \
       "-mcmc " \
       "-annotations " + ",".join(anno)

# timing = np.array([])
timing_anno = np.array([])
for i in range(250, 300):

    np.savetxt(index + "/" + gene + "/data/input.files", np.array(['C' + str(i + 1)]),
               fmt='% s', delimiter=' ', newline='\n', header='')

    data_z = pd.read_csv(index + "/" + gene + "/data" + "/" + "C" + str(i + 1) + ".txt", header=None, sep="\t")
    np.savetxt(index + "/" + gene + "/data" + "/" + "C" + str(i + 1), data_z.values[:, 1],
               fmt='%.6f', delimiter=' ', newline='\n', header='Zscore', comments='')

    os.system("cp " + index + "/" + gene + "/data" + "/LD.txt " + index + "/" + gene + "/data" + "/" + "C" + str(
        i + 1) + ".ld")

    data_anno = pd.read_csv(index + "/" + gene + '/data/C' + str(i + 1) + '_anno.txt', sep='\t')
    data_anno.rename(columns={'SNP': 'Coding'}, inplace=True)
    data_anno['Coding'] = 1
    data_anno.to_csv(index + "/" + gene + "/data" + "/" + "C" + str(i + 1) + ".annotations", index=False, sep=" ")

    # time_start = time.time()
    # os.system(cmd1)
    # time_end = time.time()
    # timing = np.append(timing, time_end - time_start)

    time_start = time.time()
    os.system(cmd2)
    time_end = time.time()
    timing_anno = np.append(timing_anno, time_end - time_start)

    os.system("rm " + index + "/" + gene + "/data" + "/input.files")

for i in range(0, 50):

    os.system("rm " + index + "/" + gene + "/data" + "/" + "C" + str(i + 1))
    os.system("rm " + index + "/" + gene + "/data" + "/" + "C" + str(i + 1) + ".ld")
    os.system("rm " + index + "/" + gene + "/data" + "/" + "C" + str(i + 1) + ".annotations")

os.system("rm " + index + "/" + gene + "/data/LD.txt")

# np.savetxt(index + "/" + gene + '/PAINTOR' + label + '.time', timing, fmt='%4f')
np.savetxt(index + "/" + gene + '/PAINTOR_anno' + label + '.time', timing_anno, fmt='%4f')
