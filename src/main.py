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

import numpy as np

def generate_total_pcoa(biom_file, tree_file, metadata_file):
	total_matrix = pairwise1.Total_Pairwise(biom_file, tree_file)
	pcoa.PCoA_total(, biom_file, metadata_file, False)
	pcoa.PCoA_group()

if __name__ == '__main__':
	generate_pcoa('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt')