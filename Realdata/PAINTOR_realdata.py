import os
import numpy as np
import pandas as pd


Chr = '21'

ld_list = np.loadtxt('chr'+Chr+'/LD/LD_list_TG.txt', dtype=int)
ld_list = ld_list[:]


def get_ldname(chr, pos):

    return 'chr'+str(chr)+'_'+str(int(pos*1e6+1))+'_'+str(int((pos+3)*1e6+1))


for Pos in ld_list:

    anno = pd.read_csv('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '.annotations', sep=' ').columns.tolist()

    cmd1 = "PAINTOR_V3.0/PAINTOR " \
           "-input " + 'chr' + Chr + '/data/input.files ' \
           "-in " + 'chr' + Chr + "/data/ " \
           "-out " + 'chr' + Chr + "/PAINTOR_result/ " \
           "-Zhead Zscore " \
           "-LDname ld " \
           "-num_samples 315133 " \
           "-MI 5 " \
           "-mcmc" \

    cmd2 = "PAINTOR_V3.0/PAINTOR " \
           "-input " + 'chr' + Chr + '/data/input.files ' \
           "-in " + 'chr' + Chr + "/data/ " \
           "-out " + 'chr' + Chr + "/PAINTOR_result_anno/ " \
           "-Zhead Zscore " \
           "-LDname ld " \
           "-num_samples 315133 " \
           "-MI 5 " \
           "-mcmc " \
           "-annotations " + ",".join(anno)

    for trait in ['HDL', 'LDL', 'Cho']:

        # trait = 'TG'

        if not os.path.exists('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '_' + trait + '.txt'):

            continue

        else:
    
            np.savetxt('chr' + Chr + '/data/input.files', np.array([get_ldname(Chr, Pos)]),
                        fmt='% s', delimiter=' ', newline='\n', header='')

            z = pd.read_csv('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '_' + trait + '.txt', header=None, sep='\t')
            np.savetxt('chr' + Chr + '/data/' + get_ldname(Chr, Pos), z.values[:, 1],
                        fmt='%.6f', delimiter=' ', newline='\n', header='Zscore', comments='')

            os.system('timeout 3h '+ cmd1)
            os.system('timeout 3h '+ cmd2)

            if os.path.exists('chr' + Chr + '/PAINTOR_result/' + get_ldname(Chr, Pos) + '.results'):
                os.system('mv chr' + Chr + '/PAINTOR_result/' + get_ldname(Chr, Pos) + '.results chr' + Chr + '/PAINTOR_result/' + get_ldname(Chr, Pos) + '_' + trait + '.results')
            if os.path.exists('chr' + Chr + '/PAINTOR_result_anno/' + get_ldname(Chr, Pos) + '.results'):
                os.system('mv chr' + Chr + '/PAINTOR_result_anno/' + get_ldname(Chr, Pos) + '.results chr' + Chr + '/PAINTOR_result_anno/' + get_ldname(Chr, Pos) + '_' + trait + '.results')  
            os.system("rm " + 'chr' + Chr + '/data/input.files')
            os.system("rm " + 'chr' + Chr + '/data/' + get_ldname(Chr, Pos))


def get_files_with_extension(directory, extension):

    pattern = f'{directory}/**/*{extension}'
    files = glob.glob(pattern, recursive=True)
    return files


def compute_paintor_credible_sets(posteriors, ld_matrix, threshold=0.95):

    n_snps = posteriors.shape[0]
    credible_sets = np.zeros(n_snps, dtype=int)
    credible_set_id = 0

    while np.sum(posteriors) > threshold:

        # Find SNP with max posterior
        max_idx = np.argmax(posteriors)

        # Get indices of SNPs in LD above threshold
        ld_indexes = np.where(ld_matrix[max_idx] > 0.5)[0]

        if ld_indexes.size != 0 and np.sum(posteriors[ld_indexes]) > threshold:
            # Subset posteriors and LD
            sub_posteriors = posteriors[ld_indexes]
            sub_posteriors /= sub_posteriors.sum()

            # Find SNPs until cumulative sum reaches threshold
            credible_indices = np.argsort(sub_posteriors)[
                               -(np.where(np.cumsum(np.sort(sub_posteriors)[::-1]) < threshold)[0].shape[0] + 1):]

            # cumsum = np.cumsum(sub_posteriors)
            # credible_indices = np.where(cumsum >= threshold)[0][0]

            credible_set_id += 1
            credible_sets[ld_indexes[credible_indices]] = credible_set_id
            # credible_sets[ld_indexes[credible_indices]] = 1

        # Remove subsetted SNPs
        # posteriors = np.delete(posteriors, ld_indexes)
        # ld_matrix = np.delete(ld_matrix, ld_indexes, axis=0)
        # ld_matrix = np.delete(ld_matrix, ld_indexes, axis=1)
        posteriors[ld_indexes] = 0

    return credible_sets


pip_load = get_files_with_extension('chr' + Chr + '/PAINTOR_result_anno/', 'TG.results')

for load in pip_load:

    LD_load = load.replace('PAINTOR_result_anno',
                            'data').replace('_HDL.results',
                                            '.ld').replace('_LDL.results',
                                                            '.ld').replace('_Cho.results', '.ld').replace('_TG.results', '.ld')
    LD = pd.read_csv(LD_load, header=None, sep=' ').values

    paintor_pip = pd.read_csv(load, sep=' ').values[:, -1]
    sets_load = load.replace('PAINTOR_result_anno', 'PAINTOR_anno_sets').replace('.results', '.sets')
    np.savetxt(sets_load, compute_paintor_credible_sets(paintor_pip, LD), fmt='%d')

    paintor_pip = pd.read_csv(load.replace('PAINTOR_result_anno', 'PAINTOR_result'), sep=' ').values[:, -1]
    sets_load = load.replace('PAINTOR_result_anno', 'PAINTOR_sets').replace('.results', '.sets')
    np.savetxt(sets_load, compute_paintor_credible_sets(paintor_pip, LD), fmt='%d')
