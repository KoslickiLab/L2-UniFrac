import sys, os
sys.path.append('../')
sys.path.append('../src')
sys.path.append('../scripts')
from os import path
import numpy as np
import matplotlib.pyplot as plt
import argparse
from argparse import ArgumentTypeError
import time
import L2Unifrac as L2U
import multiprocessing as mp
import PCoA_analysis as pcoa
import L1_pairwise_unifrac as pairwise1
import L2_pairwise_unifrac as pairwise2
import MetadataWrapper as meta
import averages as avg
import clustering as cluster
import CSVWrapper as CSV
import preprocessing as prep
import kronaUtils as krona
import diffabUtils as diff

def generate_total_pcoa(biom_file, tree_file, metadata_file, verbose, threads, intermediate_store, preprocessed_use, unifrac_code, output_file):
	if unifrac_code == 1 or unifrac_code == 2:
		if preprocessed_use and path.exists('intermediate/L1_distance_matrix_intermediate.txt'):
			total_matrix_L1 = CSV.read('intermediate/L1_distance_matrix_intermediate.txt')
			if verbose:
				print('\tSuccessfully retrieved intermediate file for L1 Generate Total PCoA')
		else:
			if verbose and preprocessed_use:
				print('\tWarning: Intermediate selected but not available. Computing pairwise alignments... This may take a while...')
			elif verbose:
				print('\tWarning: Intermediate pairwise alignment computation starting... This may take a while...')
			total_matrix_L1 = pairwise1.Total_Pairwise(biom_file, tree_file)
			if verbose:
				print('\tCompleted pairwise distance matrix computation')
			if intermediate_store:
				if verbose:
					print('\tStoring pairwise alignments...')
				if path.exists('intermediate/L1_distance_matrix_intermediate.txt'):
					os.remove('intermediate/L1_distance_matrix_intermediate.txt')
				for i in range(len(total_matrix_L1)):
					CSV.write('intermediate/L1_distance_matrix_intermediate.txt', total_matrix_L1[i])
				if verbose:
					print('\tL1 pairwise distance matrix stored successfully')
		if verbose:
			print('\tGenerating PCoA...')
		pcoa_out_L1 = pcoa.PCoA_total_from_matrix(total_matrix_L1, biom_file, metadata_file)
		if verbose:
			print('\tGeneration complete. Saving...')
		plt.savefig('images/L1_' + str(output_file) + '.png')
		if verbose:
			print('\tL1 PCoA successfully saved')
	if unifrac_code == 0 or unifrac_code == 1:
		if preprocessed_use and path.exists('intermediate/L2_distance_matrix_intermediate.txt'):
			total_matrix_L2 = CSV.read('intermediate/L2_distance_matrix_intermediate.txt')
			if verbose:
				print('\tSuccessfully retrieved intermediate file for L2 Generate Total PCoA')
		else:
			if verbose and preprocessed_use:
				print('\tWarning: Intermediate selected but not available. Computing pairwise alignments... This may take a while...')
			elif verbose:
				print('\tWarning: Intermediate pairwise alignment computation starting... This may take a while...')
			total_matrix_L2 = pairwise2.Total_Pairwise(biom_file, tree_file)
			if verbose:
				print('\tCompleted pairwise distance matrix computation')
			if intermediate_store:
				if verbose:
					print('\tStoring pairwise alignments...')
				if path.exists('intermediate/L2_distance_matrix_intermediate.txt'):
					os.remove('intermediate/L2_distance_matrix_intermediate.txt')
				for i in range(len(total_matrix_L2)):
					CSV.write('intermediate/L2_distance_matrix_intermediate.txt', total_matrix_L2[i])
				if verbose:
					print('\tL2 pairwise distance matrix stored successfully')
		if verbose:
			print('\tGenerating PCoA...')
		pcoa_out_L2 = pcoa.PCoA_total_from_matrix(total_matrix_L2, biom_file, metadata_file)
		if verbose:
			print('\tGeneration complete. Saving...')
		plt.savefig('images/L2_' + str(output_file) + '.png')
		if verbose:
			print('\tL2 PCoA successfully saved')

def generate_group_pcoa(biom_file, tree_file, metadata_file, tax_file, verbose, threads, intermediate_store, preprocessed_use, unifrac_code, output_file):
	if verbose:
		print('\tExtracting metadata...')
	metadata = meta.extract_metadata(metadata_file)
	sample_groups = []
	groups_temp = list(metadata.values())
	groups = []
	for i in range(len(groups_temp)):
		if groups_temp[i]['body_site'] not in groups:
			groups.append(groups_temp[i]['body_site'])
	group_str = ','.join(groups)
	if verbose:
		print('\tSuccessfully extracted metadata')
	if preprocessed_use and path.exists('intermediate/L1_preprocessed_intermediate.txt') and path.exists('intermediate/L2_preprocessed_intermediate.txt'):
		L1_preprocessed = CSV.read('intermediate/L1_preprocessed_intermediate.txt')
		L2_preprocessed = CSV.read('intermediate/L2_preprocessed_intermediate.txt')
		if verbose:
			print('\tSuccessfully retrieved intermediate file for L1 and L2 Preprocessing')
	else:
		if verbose and preprocessed_use:
			print('\tWarning: Intermediate selected but not available. Starting preprocessing... This may take a while...')
		elif verbose:
			print('\tWarning: Biom preprocessing starting... This may take a while...')
		if intermediate_store:
			if path.exists('intermediate/L1_preprocessed_intermediate.txt'):
				os.remove('intermediate/L1_preprocessed_intermediate.txt')
			if path.exists('intermediate/L2_preprocessed_intermediate.txt'):
				os.remove('intermediate/L2_preprocessed_intermediate.txt')
			L1_preprocessed, L2_preprocessed = prep.generate_preprocessed(biom_file, tree_file, 1, 'intermediate/L1_preprocessed_intermediate.txt', 'intermediate/L2_preprocessed_intermediate.txt')
		else:
			L1_preprocessed, L2_preprocessed = prep.generate_preprocessed(biom_file, tree_file, unifrac_code, 'tmp_L1_preprocessed_intermediate.txt', 'tmp_L2_preprocessed_intermediate.txt')
		if verbose:
			print('\tCompleted biom preprocessing matrix computation')
	if unifrac_code == 1 or unifrac_code == 2:
		if preprocessed_use or intermediate_store:
			L1_region_names, L1_tax_arr, L1_group_averages, L1_inverse_pushed, L1_neg_arr, L1_distance_matrix, L1_node_type_group_abundances = avg.compute_L1_averages('intermediate/L1_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv')
		else:
			L1_region_names, L1_tax_arr, L1_group_averages, L1_inverse_pushed, L1_neg_arr, L1_distance_matrix, L1_node_type_group_abundances = avg.compute_L1_averages('tmp_L1_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv')
		pcoa_out_L1 = pcoa.PCoA_group_from_matrix(L1_distance_matrix, biom_file, groups, plot=False)
		plt.savefig('images/out_L1_group_average.png')
	if unifrac_code == 0 or unifrac_code == 1:
		if preprocessed_use or intermediate_store:
			L2_region_names, L2_tax_arr, L2_group_averages, L2_inverse_pushed, L2_neg_arr, L2_distance_matrix, L2_node_type_group_abundances = avg.compute_L2_averages('intermediate/L2_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv')
		else:
			L2_region_names, L2_tax_arr, L2_group_averages, L2_inverse_pushed, L2_neg_arr, L2_distance_matrix, L2_node_type_group_abundances = avg.compute_L2_averages('tmp_L2_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv')
		pcoa_out_L2 = pcoa.PCoA_group_from_matrix(L2_distance_matrix, biom_file, groups, plot=False)
		plt.savefig('images/out_L2_group_average.png')
	if path.exists('tmp_L1_preprocessed_intermediate.txt'):
		os.remove('tmp_L1_preprocessed_intermediate.txt')
	if path.exists('tmp_L2_preprocessed_intermediate.txt'):
		os.remove('tmp_L2_preprocessed_intermediate.txt')

def generate_clustering_report(biom_file, tree_file, metadata_file, verbose, threads, intermediate_store, preprocessed_use, unifrac_code, output_file):
	num_clusters = meta.extract_num_clusters(metadata_file)
	if unifrac_code == 1 or unifrac_code == 2:
		if preprocessed_use and path.exists('intermediate/L1_distance_matrix_intermediate.txt'):
			total_matrix_L1 = CSV.read('intermediate/L1_distance_matrix_intermediate.txt')
			if verbose:
				print('\tSuccessfully retrieved intermediate file for L1 Generate Total PCoA')
		else:
			if verbose and preprocessed_use:
				print('\tWarning: Intermediate selected but not available. Computing pairwise alignments... This may take a while...')
			elif verbose:
				print('\tWarning: Intermediate pairwise alignment computation starting... This may take a while...')
			total_matrix_L1 = pairwise1.Total_Pairwise(biom_file, tree_file)
			if verbose:
				print('\tCompleted pairwise distance matrix computation')
			if intermediate_store:
				if verbose:
					print('\tStoring pairwise alignments...')
				if path.exists('intermediate/L1_distance_matrix_intermediate.txt'):
					os.remove('intermediate/L1_distance_matrix_intermediate.txt')
				for i in range(len(total_matrix_L1)):
					CSV.write('intermediate/L1_distance_matrix_intermediate.txt', total_matrix_L1[i])
				if verbose:
					print('\tL1 pairwise distance matrix stored successfully')
		if verbose:
			print('\tGenerating Clustering Report...')
		report = cluster.report_clustering(total_matrix_L1, biom_file, metadata_file, num_clusters, False, 1, 'reports/L1_' + str(output_file) + '_clustering.txt')
		if verbose:
			print('\tGeneration complete. Saving...')
		for i in range(len(report)):
			CSV.write('reports/L1_' + str(output_file) + '_clustering.csv', report[i])
		if verbose:
			print('\tL1 clustering successfully saved')
			print('\tGenerating Clustering PCoA...')
		pcoa = cluster.generate_clustering_pcoa(total_matrix_L1, biom_file, metadata_file, num_clusters, output_file, False, 1)
		if verbose:
			print('\tGeneration complete and saved')
	if unifrac_code == 0 or unifrac_code == 1:
		if preprocessed_use and path.exists('intermediate/L2_distance_matrix_intermediate.txt'):
			total_matrix_L2 = CSV.read('intermediate/L2_distance_matrix_intermediate.txt')
			if verbose:
				print('\tSuccessfully retrieved intermediate file for L2 Generate Total PCoA')
		else:
			if verbose and preprocessed_use:
				print('\tWarning: Intermediate selected but not available. Computing pairwise alignments... This may take a while...')
			elif verbose:
				print('\tWarning: Intermediate pairwise alignment computation starting... This may take a while...')
			total_matrix_L2 = pairwise2.Total_Pairwise(biom_file, tree_file)
			if verbose:
				print('\tCompleted pairwise distance matrix computation')
			if intermediate_store:
				if verbose:
					print('\tStoring pairwise alignments...')
				if path.exists('intermediate/L2_distance_matrix_intermediate.txt'):
					os.remove('intermediate/L2_distance_matrix_intermediate.txt')
				for i in range(len(total_matrix_L2)):
					CSV.write('intermediate/L2_distance_matrix_intermediate.txt', total_matrix_L2[i])
				if verbose:
					print('\tL2 pairwise distance matrix stored successfully')
		if verbose:
			print('\tGenerating Clustering Report...')
		report = cluster.report_clustering(total_matrix_L2, biom_file, metadata_file, num_clusters, False, 2, 'reports/L2_' + str(output_file) + '_clustering.txt')
		if verbose:
			print('\tGeneration complete. Saving...')
		for i in range(len(report)):
			CSV.write('reports/L2_' + str(output_file) + '_clustering.csv', report[i])
		if verbose:
			print('\tL2 clustering successfully saved')
			print('\tGenerating Clustering PCoA...')
		pcoa = cluster.generate_clustering_pcoa(total_matrix_L2, biom_file, metadata_file, num_clusters, output_file, False, 2)
		if verbose:
			print('\tGeneration complete and saved')

def generate_krona(biom_file, tree_file, metadata_file, tax_file, verbose, threads, intermediate_store, preprocessed_use, unifrac_code, output_file):
	if verbose:
		print('\tExtracting metadata...')
	metadata = meta.extract_metadata(metadata_file)
	sample_groups = []
	groups_temp = list(metadata.values())
	groups = []
	for i in range(len(groups_temp)):
		if groups_temp[i]['body_site'] not in groups:
			groups.append(groups_temp[i]['body_site'])
	group_str = ','.join(groups)
	if verbose:
		print('\tSuccessfully extracted metadata')
	if preprocessed_use and path.exists('intermediate/L1_preprocessed_intermediate.txt') and path.exists('intermediate/L2_preprocessed_intermediate.txt'):
		L1_preprocessed = CSV.read('intermediate/L1_preprocessed_intermediate.txt')
		L2_preprocessed = CSV.read('intermediate/L2_preprocessed_intermediate.txt')
		if verbose:
			print('\tSuccessfully retrieved intermediate file for L1 and L2 Preprocessing')
	else:
		if verbose and preprocessed_use:
			print('\tWarning: Intermediate selected but not available. Starting preprocessing... This may take a while...')
		elif verbose:
			print('\tWarning: Biom preprocessing starting... This may take a while...')
		if intermediate_store:
			if path.exists('intermediate/L1_preprocessed_intermediate.txt'):
				os.remove('intermediate/L1_preprocessed_intermediate.txt')
			if path.exists('intermediate/L2_preprocessed_intermediate.txt'):
				os.remove('intermediate/L2_preprocessed_intermediate.txt')
			L1_preprocessed, L2_preprocessed = prep.generate_preprocessed(biom_file, tree_file, 1, 'intermediate/L1_preprocessed_intermediate.txt', 'intermediate/L2_preprocessed_intermediate.txt')
		else:
			L1_preprocessed, L2_preprocessed = prep.generate_preprocessed(biom_file, tree_file, unifrac_code, 'tmp_L1_preprocessed_intermediate.txt', 'tmp_L2_preprocessed_intermediate.txt')
		if verbose:
			print('\tCompleted biom preprocessing matrix computation')
	if unifrac_code == 1 or unifrac_code == 2:
		if preprocessed_use or intermediate_store:
			L1_region_names, L1_tax_arr, L1_group_averages, L1_inverse_pushed, L1_neg_arr, L1_distance_matrix, L1_node_type_group_abundances = avg.compute_L1_averages('intermediate/L1_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv')
		else:
			L1_region_names, L1_tax_arr, L1_group_averages, L1_inverse_pushed, L1_neg_arr, L1_distance_matrix, L1_node_type_group_abundances = avg.compute_L1_averages('tmp_L1_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv')
		krona.generate_krona_visuals(L1_region_names, L1_tax_arr, L1_inverse_pushed, 'L1' + output_file, intermediate_store)
	if unifrac_code == 0 or unifrac_code == 1:
		if preprocessed_use or intermediate_store:
			L2_region_names, L2_tax_arr, L2_group_averages, L2_inverse_pushed, L2_neg_arr, L2_distance_matrix, L2_node_type_group_abundances = avg.compute_L2_averages('intermediate/L2_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv')
		else:
			L2_region_names, L2_tax_arr, L2_group_averages, L2_inverse_pushed, L2_neg_arr, L2_distance_matrix, L2_node_type_group_abundances = avg.compute_L2_averages('tmp_L2_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv')
		krona.generate_krona_visuals(L2_region_names, L2_tax_arr, L2_inverse_pushed, 'L2' + output_file, intermediate_store)
	if path.exists('tmp_L1_preprocessed_intermediate.txt'):
		os.remove('tmp_L1_preprocessed_intermediate.txt')
	if path.exists('tmp_L2_preprocessed_intermediate.txt'):
		os.remove('tmp_L2_preprocessed_intermediate.txt')

def generate_diffab(biom_file, tree_file, metadata_file, tax_file, verbose, threads, intermediate_store, preprocessed_use, unifrac_code, output_file):
	if verbose:
		print('\tExtracting metadata...')
	(Tint, lint, nodes_in_order) = L2U.parse_tree_file(tree_file)
	metadata = meta.extract_metadata(metadata_file)
	sample_groups = []
	groups_temp = list(metadata.values())
	groups = []
	for i in range(len(groups_temp)):
		if groups_temp[i]['body_site'] not in groups:
			groups.append(groups_temp[i]['body_site'])
	group_str = ','.join(groups)
	if verbose:
		print('\tSuccessfully extracted metadata')
	if preprocessed_use and path.exists('intermediate/L1_preprocessed_intermediate.txt') and path.exists('intermediate/L2_preprocessed_intermediate.txt'):
		L1_preprocessed = CSV.read('intermediate/L1_preprocessed_intermediate.txt')
		L2_preprocessed = CSV.read('intermediate/L2_preprocessed_intermediate.txt')
		if verbose:
			print('\tSuccessfully retrieved intermediate file for L1 and L2 Preprocessing')
	else:
		if verbose and preprocessed_use:
			print('\tWarning: Intermediate selected but not available. Starting preprocessing... This may take a while...')
		elif verbose:
			print('\tWarning: Biom preprocessing starting... This may take a while...')
		if intermediate_store:
			if path.exists('intermediate/L1_preprocessed_intermediate.txt'):
				os.remove('intermediate/L1_preprocessed_intermediate.txt')
			if path.exists('intermediate/L2_preprocessed_intermediate.txt'):
				os.remove('intermediate/L2_preprocessed_intermediate.txt')
			L1_preprocessed, L2_preprocessed = prep.generate_preprocessed(biom_file, tree_file, 1, 'intermediate/L1_preprocessed_intermediate.txt', 'intermediate/L2_preprocessed_intermediate.txt')
		else:
			L1_preprocessed, L2_preprocessed = prep.generate_preprocessed(biom_file, tree_file, unifrac_code, 'tmp_L1_preprocessed_intermediate.txt', 'tmp_L2_preprocessed_intermediate.txt')
		if verbose:
			print('\tCompleted biom preprocessing matrix computation')
	if unifrac_code == 1 or unifrac_code == 2:
		if preprocessed_use or intermediate_store:
			L1_region_names, L1_tax_arr, L1_group_averages, L1_inverse_pushed, L1_neg_arr, L1_distance_matrix, L1_node_type_group_abundances = avg.compute_L1_averages('intermediate/L1_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv', most_shared=True)
		else:
			L1_region_names, L1_tax_arr, L1_group_averages, L1_inverse_pushed, L1_neg_arr, L1_distance_matrix, L1_node_type_group_abundances = avg.compute_L1_averages('tmp_L1_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv', most_shared=True)
		diff.generate_diffab(L1_region_names, L1_inverse_pushed, Tint, lint, nodes_in_order, L1_tax_arr, 'L1' + output_file, 0.000005, 10, True, 1)
	if unifrac_code == 0 or unifrac_code == 1:
		if preprocessed_use or intermediate_store:
			L2_region_names, L2_tax_arr, L2_group_averages, L2_inverse_pushed, L2_neg_arr, L2_distance_matrix, L2_node_type_group_abundances = avg.compute_L2_averages('intermediate/L2_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv', most_shared=True)
		else:
			L2_region_names, L2_tax_arr, L2_group_averages, L2_inverse_pushed, L2_neg_arr, L2_distance_matrix, L2_node_type_group_abundances = avg.compute_L2_averages('tmp_L2_preprocessed_intermediate.txt', biom_file, tree_file, metadata_file, tax_file, 'reports/' + str(output_file) + '_avg_report.csv', most_shared=True)
		diff.generate_diffab(L2_region_names, L2_inverse_pushed, Tint, lint, nodes_in_order, L2_tax_arr, 'L2' + output_file, 0.000005, 10, True, 2)
	if path.exists('tmp_L1_preprocessed_intermediate.txt'):
		os.remove('tmp_L1_preprocessed_intermediate.txt')
	if path.exists('tmp_L2_preprocessed_intermediate.txt'):
		os.remove('tmp_L2_preprocessed_intermediate.txt')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='This script calculates various metrics using L2-UniFrac (or optionally L1 for comparison). '
					'It runs requests E2E without the user needing to perform intermediate steps. ', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-v', '--verbose', action='store_true', help='Print out progress report/timing information')
	
	parser.add_argument('-th', '--threads', type=int, metavar='', help='Number of threads to use',
						default=int(mp.cpu_count()/4))
	
	parser.add_argument('-i', '--intermediate_store', action='store_true', help='Stores intermediate report generation')
	parser.add_argument('-p', '--preprocessed_use', action='store_true', help='Uses intermediate report data to accelerate output')

	parser.add_argument('-bi', '--biom_file', metavar='', help='Biom file: Sample dataset in the BIOM file format')
	parser.add_argument('-tr', '--tree_file', metavar='', help='Tree file: OTU Tree in TREE file format')
	parser.add_argument('-md', '--metadata_file', metavar='', help='Metadata file: Sample names will be the first column and body habitats are the third column, following "OBERON:"')
	parser.add_argument('-tx', '--taxonomy_file', metavar='', help='Taxonomy file: Conversion between Tree Node and Taxonomy')
	
	parser.add_argument('-o', '--out_file', metavar='', help='Generic output file name. Specific requested operations will append tags and filetypes to the end')

	parser.add_argument('-t', '--total_pcoa', action='store_true', help='Compute total PCoA and generate plot')
	parser.add_argument('-g', '--group_pcoa', action='store_true', help='Compute group PCoA averages and generate plot')
	parser.add_argument('-c', '--clustering_report', action='store_true', help='Compute clustering report for distance matrix')
	parser.add_argument('-k', '--krona_averages', action='store_true', help='Compute Krona visualization for each averaged group')
	parser.add_argument('-d', '--differential_abundances', action='store_true', help='Compute differential abundances for each group pair')
	
	parser.add_argument('-L1', '--L1_UniFrac', action='store_true', help='Compute using L1 UniFrac only')
	parser.add_argument('-U', '--UniFrac', action='store_true', help='Compute using both L1 and L2 UniFrac')

	#generate_group_pcoa('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt', '../data/taxonomies/gg_13_8_99.gg.tax')
	#generate_total_pcoa('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt')

	start = time.time()
	args = parser.parse_args()
	print(args)

	if ((isinstance(args.biom_file, str) and not path.exists(args.biom_file)) or args.biom_file is None) and (args.total_pcoa or args.group_pcoa or args.clustering_report):
		raise Exception('Error: Invalid Biom File Path')
	elif ((isinstance(args.tree_file, str) and not path.exists(args.tree_file)) or args.tree_file is None) and (args.total_pcoa or args.group_pcoa or args.clustering_report):
		raise Exception('Error: Invalid Tree File Path')
	elif ((isinstance(args.metadata_file, str) and not path.exists(args.metadata_file)) or args.metadata_file is None) and (args.total_pcoa or args.group_pcoa or args.clustering_report):
		raise Exception('Error: Invalid Metadata File Path')
	elif ((isinstance(args.taxonomy_file, str) and not path.exists(args.taxonomy_file)) or args.taxonomy_file is None) and (args.group_pcoa):
		raise Exception('Error: Invalid Taxonomy File Path')
	elif args.out_file is None and (args.total_pcoa or args.group_pcoa or args.clustering_report):
		raise Exception('Error: Missing output file subname. Please specify and continue')

	unifrac_code = 0
	if args.UniFrac:
		unifrac_code = 1
	elif args.L1_UniFrac:
		unifrac_code = 2

	if path.exists(args.out_file) and args.verbose:
		print('Output file already exists... It will be overwritten.')

	if args.total_pcoa:
		segment_start = time.time()
		if args.verbose:
			print('Starting Total PCoA Generation...')
		generate_total_pcoa(args.biom_file, args.tree_file, args.metadata_file, args.verbose, args.threads, args.intermediate_store, args.preprocessed_use, unifrac_code, args.out_file)
		if args.verbose:
			print('Total PCoA Generation Complete. Total Elapsed Time: ' + str(time.time()-segment_start) + ' seconds')
	if args.group_pcoa:
		segment_start = time.time()
		if args.verbose:
			print('Starting Group PCoA Generation...')
		generate_group_pcoa(args.biom_file, args.tree_file, args.metadata_file, args.taxonomy_file, args.verbose, args.threads, args.intermediate_store, args.preprocessed_use, unifrac_code, args.out_file)
		if args.verbose:
			print('Group PCoA Generation Complete. Total Elapsed Time: ' + str(time.time()-segment_start) + ' seconds')
	if args.clustering_report:
		segment_start = time.time()
		if args.verbose:
			print('Starting Clustering Report Generation...')
		generate_clustering_report(args.biom_file, args.tree_file, args.metadata_file, args.verbose, args.threads, args.intermediate_store, args.preprocessed_use, unifrac_code, args.out_file)
		if args.verbose:
			print('Clustering Report Generation Complete. Total Elapsed Time: ' + str(time.time()-segment_start) + ' seconds')
	if args.krona_averages:
		segment_start = time.time()
		if args.verbose:
			print('Starting Krona Visualization Generation...')
		generate_krona(args.biom_file, args.tree_file, args.metadata_file, args.taxonomy_file, args.verbose, args.threads, args.intermediate_store, args.preprocessed_use, unifrac_code, args.out_file)
		if args.verbose:
			print('Krona Visualization Generation Complete. Total Elapsed Time: ' + str(time.time()-segment_start) + ' seconds')
	if args.differential_abundances:
		segment_start = time.time()
		if args.verbose:
			print('Starting Differential Abundance Visualization Generation...')
		generate_diffab(args.biom_file, args.tree_file, args.metadata_file, args.taxonomy_file, args.verbose, args.threads, args.intermediate_store, args.preprocessed_use, unifrac_code, args.out_file)
		if args.verbose:
			print('Differential Abundance Visualization Generation Complete. Total Elapsed Time: ' + str(time.time()-segment_start) + ' seconds')
