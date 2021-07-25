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
import MetadataWrapper as meta
import TaxWrapper as tax

import PCoA_analysis as pcoa
import pairwise_unifrac as pairwise1
import L1_pairwise_unifrac as pairwise2
import matplotlib.pyplot as plt

import numpy as np

def generate_total_pcoa(biom_file, tree_file, metadata_file):
	total_matrix_L1 = pairwise2.Total_Pairwise(biom_file, tree_file)
	total_matrix_L2 = pairwise1.Total_Pairwise(biom_file, tree_file)
	pcoa_out_L1 = pcoa.PCoA_total_from_matrix(total_matrix_L1, biom_file, metadata_file, False)
	plt.savefig('images/out_L1.png')
	pcoa_out_L1 = pcoa.PCoA_total_from_matrix(total_matrix_L2, biom_file, metadata_file, False)
	plt.savefig('images/out_L2.png')

if __name__ == '__main__':
	generate_total_pcoa('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt')