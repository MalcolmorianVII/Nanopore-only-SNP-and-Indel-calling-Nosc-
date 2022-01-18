import pandas as pd
import os.path
import csv

clair = '/Users/malcolmorian/Documents/Bioinformatics/Projects2021/Guppy3Guppy5/NOSC/nosc_clair/2022.01.02/clair_vcfData'
pepper = '/Users/malcolmorian/Documents/Bioinformatics/Projects2021/Guppy3Guppy5/NOSC/nosc_pepper/2022.01.02/pepper_vcfData'
gatk = '/Users/malcolmorian/Documents/Bioinformatics/Projects2021/Guppy3Guppy5/NOSC/nosc_gatk/2022.01.05/gatk_vcfData'

venned ='/Users/malcolmorian/Documents/Bioinformatics/Projects2021/Guppy3Guppy5/NOSC/vennTextFiles'
# os.mkdir(venned)


def write(data_df,toWrite):
    with open(toWrite,'w') as f:
        # print(toWrite)
        # if os.path.exists(toWrite):
        #     print('NOT WRITING ....THE FILE ALREADY EXISTS')
        data_df.to_csv(toWrite,sep=' ', index=False, header=False)


def read_file(directory):
    files = os.listdir(directory)
    print(files)
    for file in files:
        if file.endswith('.xlsx'):
            path = f'{directory}/{file}'
            print(path)
            data_df = pd.read_excel(path, engine='openpyxl')
            data_df = data_df[['#CHROM','POS']]

            write(data_df, f'{venned}/{file}.csv')
        else:
            print('NOT THIS ONE!!!!!!!!!!!')
            continue


read_file(clair)
read_file(gatk)
read_file(pepper)

