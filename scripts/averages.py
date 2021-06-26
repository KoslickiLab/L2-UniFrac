import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import CSVWrapper as CSV
import MetadataWrapper as meta
import numpy as np

# Note: these are the same for L1/L2, so they will be computed only once.
nodes_samples = BW.extract_biom('../data/47422_otu_table.biom')
T1, l1, nodes_in_order = L2U.parse_tree_file('../data/trees/gg_13_5_otus_99_annotated.tree')
(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')
metadata = meta.extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt')

region_names = []
region_map = {}
for i in range(len(PCoA_Samples)):
	if metadata[PCoA_Samples[i]]['body_site'] not in region_names:
		region_map[metadata[PCoA_Samples[i]]['body_site']] = []
		region_names.append(metadata[PCoA_Samples[i]]['body_site'])
	region_map[metadata[PCoA_Samples[i]]['body_site']].append(i)
	PCoA_Samples[i] = region_names.index(metadata[PCoA_Samples[i]]['body_site'])

sparse_matrix_L1 = CSV.read_sparse('L1-Push-Out.csv')
sparse_matrix_L2 = CSV.read_sparse('L2-Push-Out.csv')

group_averages_L1 = {}
group_averages_L2 = {}

CSV.write('Group-Averages.csv', region_names)

for i in range(len(region_names)):
	group_arr = []
	for j in range(len(region_map[region_names[i]])):
		group_arr.append(np.array(sparse_matrix_L1[region_map[region_names[i]][j]].todense())[0])
	average = L1U.median_of_vectors(group_arr)
	group_averages_L1[region_names[i]] = average

print("L1 Group Averages:")
CSV.write('Group-Averages.csv', ["L1 Group Averages:"])
for name in region_names:
	padded_name = "{:<15}".format(name+":")
	print(f"{padded_name} {group_averages_L1[name]}")
	CSV.write('Group-Averages.csv', group_averages_L1[name])

for i in range(len(region_names)):
	group_arr = []
	for j in range(len(region_map[region_names[i]])):
		group_arr.append(np.array(sparse_matrix_L2[region_map[region_names[i]][j]].todense())[0])
	average = L2U.mean_of_vectors(group_arr)
	group_averages_L2[region_names[i]] = average

print("\nL2 Group Averages:")
CSV.write('Group-Averages.csv', ["L2 Group Averages:"])
for name in region_names:
	padded_name = "{:<15}".format(name+":")
	print(f"{padded_name} {group_averages_L2[name]}")
	CSV.write('Group-Averages.csv', group_averages_L2[name])

print("\nL1 Pushed Up:")
CSV.write('Group-Averages.csv', ["L1 Pushed Up:"])
L1_neg_arr = []
for name in region_names:
	neg_count = 0
	median_inverse = L1U.inverse_push_up(group_averages_L1[name], T1, l1, nodes_in_order)
	for i in range(len(median_inverse)):
		if median_inverse[i] < 0:
			neg_count += 1
	L1_neg_arr.append(neg_count)
	padded_name = "{:<15}".format(name+":")
	print(f"{padded_name} {median_inverse}")
	CSV.write('Group-Averages.csv', median_inverse)

print("\nL2 Pushed Up:")
CSV.write('Group-Averages.csv', ["L2 Pushed Up:"])
L2_neg_arr = []
for name in region_names:
	neg_count = 0
	mean_inverse = L2U.inverse_push_up(group_averages_L2[name], T1, l1, nodes_in_order)
	for i in range(len(mean_inverse)):
		if mean_inverse[i] < 0:
			neg_count += 1
	L2_neg_arr.append(neg_count)
	padded_name = "{:<15}".format(name+":")
	print(f"{padded_name} {mean_inverse}")
	CSV.write('Group-Averages.csv', mean_inverse)

CSV.write('Group-Averages.csv', ["L1 and L2 Negatives by Group:"])
CSV.write('Group-Averages.csv', L1_neg_arr)
CSV.write('Group-Averages.csv', L2_neg_arr)
