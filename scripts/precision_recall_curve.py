import itertools
import os
import csv
import pandas as pd
import glob


# Input : pep & gatk vcfs i.e. originated with illumina reads


def read_file(file):
    return pd.read_excel(file, engine='openpyxl')


def add_config():
    config_file = '/Users/malcolmorian/Documents/Bioinformatics/Projects2021/Guppy3Guppy5/NOSC/nosc_genome_length_config.txt'  # path of config file
    configs = {}
    with open(config_file, 'r') as configuration:
        config = configuration.readlines()
        for line in config:
            sp, gen_len = line.strip().split(',')
            configs.update({sp: gen_len})

    return configs


def add_rpqa(df):  # re write this
    df['RefPosAltQual'] = ""
    for index, row in df.iterrows():  # Adds a tuple of (ref_pos_alt,qual) to df
        rfpsa = row['REF'] + str(row['POS']) + row['ALT'], row['QUAL']
        df['RefPosAltQual'][index] = tuple(rfpsa, )
    return df


def generate_vcf_chunks(genome_len: int, chunk_length: int):
    """

    :param chunk_total: How many chunks per split data set
    :param genome_len: Total length of a given genome
    :param chunk_length: Specifies how to split the variant sites into chunks
    :return train_set: Contains the first half of the variant sites
    """
    genome_len = int(genome_len)
    chunk_total = round(genome_len / 100000)
    half = int(genome_len / 2)
    train_pos = set(range(1, (half + 1)))
    # test_pos = set(range((half+1),gen_len+1))
    i = 0  # Initializes the chunk position considering the genome length
    train_set = []
    for q in range(chunk_total):
        chunk = tuple(itertools.islice(train_pos, i, i + chunk_length), )
        train_set.append(chunk)
        i += chunk_length
    return train_set


def get_roc_values(tool_var, gen_len, gatk_vars):
    output = {}
    # print('tool_variants',tool_var)
    for quality in range(100):
        positives = [x for x in tool_var if x[1] >= quality]
        neg = [x for x in tool_var if x[1] < quality]  # this isn't all your negatives
        tp = len([x for x in positives if x[0] in gatk_vars])
        fp = len([x for x in positives if x[0] not in gatk_vars])
        print("fps:",fp)
        print("tps:",tp)
        all_negatives = int(gen_len) - len(positives)
        # all_negatives = neg + fp
        fn = len(gatk_vars) - tp
        tn = all_negatives - fn
        try:
            tpr = tp / (tp + fn)  # 70/(70 + 30) aka recall,sensitivity
            ppv = tp / (tp + fp)  # aka precision
            # print("TPR", tpr)
            fpr = fp / (fp + tn)
            # print("FPR", fpr)
        except ZeroDivisionError as ze:
            ppv = 0
        output[quality] = (fpr, tpr, ppv)
    return output


def write_roc_values(out_dict, outfile):
    with open(outfile, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Threshold', 'fpr', 'tpr', 'precision'])
        for q in out_dict:
            writer.writerow([q, out_dict[q][0], out_dict[q][1], out_dict[q][2]])


def get_tool_vars(rpqed_df):
    # out_df = rpqed_df['REF_POS_ALT']
    tool_vars = [x for x in rpqed_df['RefPosAltQual']]
    return tool_vars


def generate_roc_values(train_set, tool_vars, gatk_vars, genome_len, outfile):
    """
    (i) Generates pepper variants based on the train set
    (ii) Use get_roc_values to generate roc_values from pepper_variant & gatk_variant
    i.e return {'Threshold':('fpr','tpr')}
    :param train_set:
    :return:
    """
    # if os.path.exists(outfile):
    #     return 'EXIT STATUS::Already run'
    tool_vars_train = [x for x in tool_vars if int(x[0].rstrip('ATCGN-,').lstrip('ATCGN-,')) in train_set]
    # print('tools_var_train length',len(tool_vars_train))
    roc_values = get_roc_values(tool_vars_train, genome_len, gatk_vars)
    write_roc_values(roc_values, outfile)
    print('Finished generating roc values')


def runner(base_dir):
    conf = add_config()
    species = os.listdir(base_dir)
    count = 0
    for sp in species:
        count += 1
        print('count',count)
        print('Processing',sp)
        if sp == '.DS_Store':
            print('SKIPPING THIS ONE!!!!!')
            continue
        # out file
        outfile_pep = f"{base_dir}/{sp}/{sp}_prc_pepper.csv"
        outfile_clair = f"{base_dir}/{sp}/{sp}_prc_clair.csv"

        if os.path.exists(outfile_pep) and os.path.exists(outfile_clair):
            print('EXIT STATUS::Already run')
            continue

        pep_vcf = ''.join(glob.glob(f'{base_dir}/{sp}/*pepper.xlsx'))
        clair_vcf = ''.join(glob.glob(f'{base_dir}/{sp}/*clair.xlsx'))
        gatk_vcf = ''.join(glob.glob(f'{base_dir}/{sp}/*gatk.xlsx'))

        # read vcfs files &  add_rpqa
        pep_df = read_file(pep_vcf)
        clair_df = read_file(clair_vcf)
        gatk_df = read_file(gatk_vcf)

        # Add rpqa
        pep_df = add_rpqa(pep_df)
        clair_df = add_rpqa(clair_df)
        gatk_df = add_rpqa(gatk_df)

        # gatk_set
        sp_genome = conf.get(sp)
        # Generate test set i,e. call generate_vcf_chunks
        train_set = generate_vcf_chunks(sp_genome, 100000)
        train_set = set(value for x in train_set for value in x)

        # Get tool variants
        pep_vars = get_tool_vars(pep_df)
        clair_vars = get_tool_vars(clair_df)
        gatk_vars = get_tool_vars(gatk_df)
        gatk_vars = [x[0] for x in gatk_vars ]


        # Generate roc_values for both clair & pep
        generate_roc_values(train_set, pep_vars, gatk_vars, sp_genome, outfile_pep)
        generate_roc_values(train_set, clair_vars, gatk_vars, sp_genome, outfile_clair)
        # print('count:::',count)

    return '<<< Message: Run session complete >>>'


base = '/Users/malcolmorian/Documents/Bioinformatics/Projects2021/Guppy3Guppy5/NOSC/all_species'
if __name__ == "__main__":
    runner(base)
