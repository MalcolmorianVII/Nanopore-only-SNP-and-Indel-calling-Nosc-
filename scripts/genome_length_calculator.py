#! /home/belson/anaconda3/bin/python
import re
import os
import glob

def run_calculate_alternative(data_list):
	length = 0
	for index,line in enumerate(data_list):
		if line.__contains__('chromosome'):
			data_list = data_list[index+1:]
			break
	for line in data_list:
			if line.__contains__('plasmid'):
				print('Reached a plasmid! Stopping genome length calculation')
				break
			length += len(line.strip())
	return length


def read_file(file_path):
	with open(file_path,'r') as f:
		file = f.readlines()
		return file
def write_genome_length(gen_len,species,output):
	with open(output,'a') as out_file:
		out_file.write(f'{species},{gen_len}\n')

def calc_gen_length(file):
	file_data = read_file(file)
	print(f'CALCULATING THE GENOME LENGTH of {file}')
	genome_length = run_calculate_alternative(file_data)
	print('GENOME CALCULATION DONE')
	return genome_length


def runner(path):
	files = glob.glob(f'{path}/*.fna')
	outfile = f'{path}/genome_lengths.txt'
	if os.path.exists(outfile):
		os.remove(outfile)
	for file in files:
		print(f'Reading the contents of {file}')
		data = read_file(file)
		length = run_calculate_alternative(data)
		species = species = re.compile(r'[A-Z]\w+(?!_ref.fna)?').search(file).group()
		print(f'Writing genome length of {file}')
		write_genome_length(length,species,outfile)

ref='/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/reference_genomes'
if __name__ == '__main__':
	runner(ref)
	print('Finished writing the genome lengths of all species')


