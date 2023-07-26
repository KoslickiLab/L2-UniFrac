import os
import sys
sys.path.append('L2-UniFrac')
sys.path.append('L2-UniFrac/src')
sys.path.append('L2-UniFrac/scripts')
sys.path.append('src')
import argparse
import L2UniFrac as L2U
import numpy as np
import pandas as pd
from helper import get_metadata_dict, get_pheno_sample_dict, get_rep_sample_dict_wgs, write_rep_samples_to_file


def generate_rep_sample_from_metadata(meta_dict, profile_dir, outfile_name, out_format='cami', leaves_only=False):
	'''

	:param meta_dict:
	:param profile_list:
	:param save_dir:
	:return:
	'''
	profile_list = os.listdir(profile_dir)
	for profile_name in profile_list:
		if not profile_name.endswith('.profile'):
			profile_list.remove(profile_name)
	profile_name_list = list(map(lambda x: x.split('.')[0], profile_list))
	profile_path_lst = [os.path.join(profile_dir, file) for file in profile_list]
	Tint, lint, nodes_in_order, nodes_to_index = L2U.get_wgs_tree(profile_path_lst)
	targets = [meta_dict[i] for i in profile_name_list]
	#profile list should come from meta_dict
	pheno_sample_dict = get_pheno_sample_dict(profile_path_lst, targets)
	rep_sample_dict = get_rep_sample_dict_wgs(pheno_sample_dict, Tint, lint, nodes_in_order, nodes_to_index)
	if out_format == 'otu':
		write_rep_samples_to_file(rep_sample_dict, outfile_name, nodes_in_order, nodes_to_index)
	else:
		index_to_nodes = {v: k for k, v in nodes_to_index.items()}
		rep_profiles_dict = L2U.build_profiles_from_dict(rep_sample_dict, nodes_in_order, index_to_nodes, leaves_only)
		for id, profile in rep_profiles_dict.items():
			out_file_name = outfile_name.split('.')[0] + '_' + str(id) + '.profile'
			profile.write_CAMI_file(out_file_name)
	return


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Get representative samples from a metadata file.')
	parser.add_argument('-m', '--meta_file', type=str, help='A metadata file.', nargs='?')
	parser.add_argument('-id_col', '--id_col', type=str, help='Name of the id column in the metadata file.', nargs='?', default="library_id")
	parser.add_argument('-s', '--save', type=str, help="Save output file as.")
	parser.add_argument('-d', '--pdir', type=str, help="Directory of profiles")
	parser.add_argument('-l', '--leaves_only', type=str, help="Only have abundances on leaf notes?", nargs='?', default='n', choices=['y','n'])
	parser.add_argument('-e', '--env_col', type=str, help='A selected phenotype corresponding to a column name in the metadata file.', nargs='?', default="HMgDB_diagnosis")
	parser.add_argument('-f', '--out_format', type=str, help='The format of output files. Choices: cami, otu. If otu is chosen, '
														'-s flag is required. Otherwise, -o is required.', choices=['cami', 'otu'], nargs='?', default='cami')

	args = parser.parse_args()
	metadata_file = args.meta_file
	profile_dir = args.pdir
	metadata_key = args.phenotype
	id_col = args.id_col
	outfile_name = args.save
	if args.leaves_only == 'y':
		leaves_only = True
	else:
		leaves_only = False

	meta_dict = get_metadata_dict(metadata_file, val_col=metadata_key, key_col=id_col)
	rep_sample_dict = generate_rep_sample_from_metadata(meta_dict, profile_dir, outfile_name, args.out_format, args.leaves_only)


