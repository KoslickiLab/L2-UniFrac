# The purpose of main.py is to perform all necessary operations using a base set of data.

import sys
sys.path.append('../')
sys.path.append('../src')
sys.path.append('../scripts')
from os import path
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import CSVWrapper as CSV
import TaxWrapper as tax

import PCoA_analysis as pcoa
import L1_pairwise_unifrac as pairwise1
import pairwise_unifrac as pairwise2
import matplotlib.pyplot as plt
import MetadataWrapper as meta
import averages as avg
import clustering as cluster
import argparse
from argparse import ArgumentTypeError
import multiprocessing as mp

import numpy as np

def generate_total_pcoa(biom_file, tree_file, metadata_file, verbose=False, intermediate_store=False, preprocessed_use=False):
	total_matrix_L1 = pairwise1.Total_Pairwise(biom_file, tree_file)
	total_matrix_L2 = pairwise2.Total_Pairwise(biom_file, tree_file)
	pcoa_out_L1 = pcoa.PCoA_total_from_matrix(total_matrix_L1, biom_file, metadata_file)
	plt.savefig('images/out_L1.png')
	pcoa_out_L1 = pcoa.PCoA_total_from_matrix(total_matrix_L2, biom_file, metadata_file)
	plt.savefig('images/out_L2.png')

def generate_group_pcoa(biom_file, tree_file, metadata_file, tax_file, verbose=False, intermediate_store=False, preprocessed_use=False):
	metadata = meta.extract_metadata(metadata_file)
	sample_groups = []
	groups_temp = list(metadata.values())
	groups = []
	for i in range(len(groups_temp)):
		if groups_temp[i]['body_site'] not in groups:
			groups.append(groups_temp[i]['body_site'])
	group_str = ','.join(groups)
	if not preprocessed_use:
		L1_preprocessed, L2_preprocessed = generate_preprocessed(biom_file, tree_file)
	_, _, _, _, _, _, _, _, L1_distance_matrix, L2_distance_matrix, _, _ = avg.compute_averages(L1_preprocessed, L2_preprocessed, biom_file, tree_file, metadata_file, tax_file)
	pcoa_out_L1 = pcoa.PCoA_group_from_matrix(L1_distance_matrix, biom_file, group_str, plot=False)
	plt.savefig('images/out_L1_group_average.png')
	pcoa_out_L2 = pcoa.PCoA_group_from_matrix(L2_distance_matrix, biom_file, group_str, plot=False)
	plt.savefig('images/out_L2_group_average.png')

def generate_clustering_report(biom_file, tree_file, metadata_file, verbose=False, intermediate_store=False, preprocessed_use=False):
	if not preprocessed_use:
		total_matrix_L1 = pairwise1.Total_Pairwise(biom_file, tree_file)
		total_matrix_L2 = pairwise2.Total_Pairwise(biom_file, tree_file)
	cluster.report_clustering(total_matrix_L1, total_matrix_L2, biom_file, metadata_file, 'reports/clustering_report.txt')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="This script calculates various metrics using L2-UniFrac (or optionally L1 for comparison). "
					"It runs requests E2E without the user needing to perform intermediate steps. ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-v', '--verbose', action="store_true", help="Print out progress report/timing information")
	
	parser.add_argument('-th', '--threads', type=int, metavar='', help="Number of threads to use",
						default=mp.cpu_count())
	
	parser.add_argument('-i', '--intermediate_store', action="store_true", help="Stores intermediate report generation")
	parser.add_argument('-p', '--preprocessed_use', action="store_true", help="Uses intermediate report data to accelerate output")

	parser.add_argument('-bi', '--biom_file', metavar='', help="Biom file: Sample dataset in the BIOM file format")
	parser.add_argument('-tr', '--tree_file', metavar='', help="Tree file: OTU Tree in TREE file format")
	parser.add_argument('-md', '--metadata_file', metavar='', help="Metadata file: Sample names will be the first column and body habitats are the third column, following 'OBERON:'")
	parser.add_argument('-tx', '--taxonomy_file', metavar='', help="Taxonomy file: Conversion between Tree Node and Taxonomy")
	
	parser.add_argument('-o', '--out_file', metavar='', help='Output csv file with the containment indices')

	parser.add_argument('-t', '--total_pcoa', action="store_true", help='Compute total PCoA and generate plot')
	parser.add_argument('-g', '--group_pcoa', action="store_true", help='Compute group PCoA averages and generate plot')
	parser.add_argument('-c', '--clustering_report', action="store_true", help='Compute clustering report for distance matrix')
	
	parser.add_argument('-L1', '--L1_UniFrac', action="store_true", help='Compute using L1 UniFrac only')
	parser.add_argument('-L2', '--L2_UniFrac', action="store_true", help='Compute using L2 UniFrac only')
	parser.add_argument('-U', '--UniFrac', action="store_true", help='Compute using both L1 and L2 UniFrac')

	#generate_group_pcoa('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt', '../data/taxonomies/gg_13_8_99.gg.tax')
	#generate_total_pcoa('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt')

	args = parser.parse_args()
	print(args.verbose)