import os
import numpy as np
import pandas as pd
from Funmap import SUSIE
from Funmap import FUNMAP

np.set_printoptions(threshold=1e8)

Chr = '21'

ld_list = np.loadtxt('chr'+Chr+'/LD/LD_list_TG.txt', dtype=int)
ld_list = ld_list[:]


def get_ldname(chr, pos):

    return 'chr'+str(chr)+'_'+str(int(pos*1e6+1))+'_'+str(int((pos+3)*1e6+1))


for Pos in ld_list:

    R = pd.read_csv('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '.ld', header=None, sep=' ')
    R = R.values

    A = pd.read_csv('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '.annotations', sep=' ')
    A = A.values.astype(float)

    n = 315133

    for trait in ['HDL', 'LDL', 'Cho']:

        # trait = 'TG'

        if not os.path.exists('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '_' + trait + '.txt'):

            continue

        else:

            z = pd.read_csv('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '_' + trait + '.txt', header=None, sep='\t')
            snp_name = z.values[:, 0]
            z = z.values[:, 1]

            result = SUSIE(z, R, n, L=10)
            s_result = pd.DataFrame({'SNP': snp_name, 'PIP': result.pip})
            s_result.to_csv('chr' + Chr + '/SuSiE_result/' + get_ldname(Chr, Pos) + '_' + trait + '.csv', header=None, index=False)
            file = open('chr' + Chr + '/SuSiE_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets', "w")
            file.write(str(result.sets))
            file.close()

            result = FUNMAP(z, R, A, n, L=10)
            s_result = pd.DataFrame({'SNP': snp_name, 'PIP': result.pip})
            s_result.to_csv('chr' + Chr + '/Funmap_result/' + get_ldname(Chr, Pos) + '_' + trait + '.csv', header=None, index=False)
            file = open('chr' + Chr + '/Funmap_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets', "w")
            file.write(str(result.sets))
            file.close()
