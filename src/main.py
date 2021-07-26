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

import numpy as np

def generate_total_pcoa(biom_file, tree_file, metadata_file):
	total_matrix_L1 = pairwise1.Total_Pairwise(biom_file, tree_file)
	total_matrix_L2 = pairwise2.Total_Pairwise(biom_file, tree_file)
	pcoa_out_L1 = pcoa.PCoA_total_from_matrix(total_matrix_L1, biom_file, metadata_file)
	plt.savefig('images/out_L1.png')
	pcoa_out_L1 = pcoa.PCoA_total_from_matrix(total_matrix_L2, biom_file, metadata_file)
	plt.savefig('images/out_L2.png')

def generate_group_pcoa(biom_file, tree_file, metadata_file, tax_file):
	metadata = meta.extract_metadata(metadata_file)
	sample_groups = []
	groups_temp = list(metadata.values())
	groups = []
	for i in range(len(groups_temp)):
		if groups_temp[i]['body_site'] not in groups:
			groups.append(groups_temp[i]['body_site'])
	group_str = ','.join(groups)
	L1_preprocessed, L2_preprocessed = generate_preprocessed(biom_file, tree_file)
	_, _, _, _, _, _, _, _, L1_distance_matrix, L2_distance_matrix, _, _ = avg.compute_averages(L1_preprocessed, L2_preprocessed, biom_file, tree_file, metadata_file, tax_file)
	pcoa_out_L1 = pcoa.PCoA_group_from_matrix(L1_distance_matrix, biom_file, group_str, plot=False)
	plt.savefig('images/out_L1_group_average.png')
	pcoa_out_L2 = pcoa.PCoA_group_from_matrix(L2_distance_matrix, biom_file, group_str, plot=False)
	plt.savefig('images/out_L2_group_average.png')

def generate_clustering_report(biom_file, tree_file, metadata_file):
	total_matrix_L1 = pairwise1.Total_Pairwise(biom_file, tree_file)
	total_matrix_L2 = pairwise2.Total_Pairwise(biom_file, tree_file)
	cluster.report_clustering(total_matrix_L1, total_matrix_L2, biom_file, metadata_file, 'reports/clustering_report.txt')

if __name__ == '__main__':
	#generate_group_pcoa('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt', '../data/taxonomies/gg_13_8_99.gg.tax')
	generate_total_pcoa('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt')