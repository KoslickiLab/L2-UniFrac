import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import write_to_csv as CSV
import MetadataWrapper as meta
import multiprocessing as mp
import write_to_csv as CSV
import numpy

cores = mp.cpu_count()

# Note: these are the same for L1/L2, so they will be computed only once.
nodes_samples = BW.extract_biom('../data/47422_otu_table.biom')
T1, l1, nodes_in_order = L2U.parse_tree_file('../data/trees/gg_13_5_otus_99_annotated.tree')
(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')

def L1_pushup_worker(sample_num):
	L1_Pushed = L1U.push_up(nodes_weighted[PCoA_Samples[sample_num]], T1, l1, nodes_in_order)
	return L1_Pushed

def L2_pushup_worker(sample_num):
	L2_Pushed = L2U.push_up(nodes_weighted[PCoA_Samples[sample_num]], T1, l1, nodes_in_order)
	return L2_Pushed

# Multi Core Method
values = range(len(PCoA_Samples))

with mp.Pool(processes=cores/2-2) as pool:
	L1_result = pool.map(L1_pushup_worker, values)

for i in range(len(L1_result)):
	for j in range(len(L1_result[i])):
		if L1_result[i][j] != 0:
			CSV.write('L1-Push-Out.csv', [i, j, L1_result[i][j]])

with mp.Pool(processes=cores/2-2) as pool:
	L2_result = pool.map(L2_pushup_worker, values)

for i in range(len(L2_result)):
	for j in range(len(L2_result[i])):
		if L2_result[i][j] != 0:
			CSV.write('L2-Push-Out.csv', [i, j, L2_result[i][j]])
