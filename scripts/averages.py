import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
from os import path
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import CSVWrapper as CSV
import MetadataWrapper as meta
import TaxWrapper as tax
import numpy as np

# File cheatsheet (python averages.py L1-Push-Out.csv L2-Push-Out.csv ../data/47422_otu_table.biom ../data/trees/gg_13_5_otus_99_annotated.tree ../data/metadata/P_1928_65684500_raw_meta.txt ../data/taxonomies/gg_13_8_99.gg.tax Group-Averages.csv):
# L1_file:       'L1-Push-Out.csv'
# L2_file:       'L2-Push-Out.csv'
# biom_file:     '../data/47422_otu_table.biom'
# tree_file:     '../data/trees/gg_13_5_otus_99_annotated.tree'
# metadata_file: '../data/metadata/P_1928_65684500_raw_meta.txt'
# tax_file:      '../data/taxonomies/gg_13_8_99.gg.tax'
# output_file:   'Group-Averages.csv'

# Negative values can periodically appear on machines when the value should be 0. This filter is
# to ignore very small negatives between -10e-14 and 0 to account for this. We tested other
# thresholds such as -10e-12 with the same results, so we know that all remaining negatives
# are very large relative to these erroneous ones.
negatives_filtering_threshold = -10e-14

def compute_pairwise_pushed(pushed_arr):
	dist_matrix = []
	for i in range(len(pushed_arr)):
		dist_arr = []
		for j in range(len(pushed_arr)):
			unifrac_distance = np.linalg.norm(pushed_arr[i] - pushed_arr[j])
			dist_arr.append(unifrac_distance)
		dist_matrix.append(dist_arr)
	return dist_matrix


# Helper function for averages. Outputs a CSV containing info from each step and negative counts at end.
def compute_averages(L1_file, L2_file, biom_file, tree_file, metadata_file, tax_file, output_file):

	# Note: these are the same for L1/L2, so they will be computed only once.
	nodes_samples = BW.extract_biom(biom_file)
	T1, l1, nodes_in_order = L2U.parse_tree_file(tree_file)
	#(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

	PCoA_Samples = BW.extract_samples(biom_file)
	metadata = meta.extract_metadata(metadata_file)

	# Extract region names and map samples to regions
	region_names = []
	region_map = {}
	for i in range(len(PCoA_Samples)):
		if metadata[PCoA_Samples[i]]['body_site'] not in region_names:
			region_map[metadata[PCoA_Samples[i]]['body_site']] = []
			region_names.append(metadata[PCoA_Samples[i]]['body_site'])
		region_map[metadata[PCoA_Samples[i]]['body_site']].append(i)
		PCoA_Samples[i] = region_names.index(metadata[PCoA_Samples[i]]['body_site'])

	# Read sparse matrices
	sparse_matrix_L1 = CSV.read_sparse(L1_file)
	sparse_matrix_L2 = CSV.read_sparse(L2_file)

	group_averages_L1 = {}
	group_averages_L2 = {}

	# Store region names for later
	print(region_names)
	CSV.write(output_file, region_names)

	# Write taxas for cell
	taxonomies = tax.extract_tax(tax_file)
	tax_arr = []
	for i in range(len(nodes_in_order)):
		if nodes_in_order[i][0] != 't':
			tax_arr.append(taxonomies[int(nodes_in_order[i])])
		else:
			tax_arr.append(nodes_in_order[i])
	print(tax_arr)
	CSV.write(output_file, tax_arr)

	# Take L1 average of each
	L1_pushed_arr = []
	for i in range(len(region_names)):
		group_arr = []
		for j in range(len(region_map[region_names[i]])):
			group_arr.append(np.array(sparse_matrix_L1[region_map[region_names[i]][j]].todense())[0])
		average = L1U.median_of_vectors(group_arr)
		group_averages_L1[region_names[i]] = average
		L1_pushed_arr.append(average)

	# Store L1 averages
	print("L1 Group Averages:")
	CSV.write(output_file, ["L1 Group Averages:"])
	for name in region_names:
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {group_averages_L1[name]}")
		CSV.write(output_file, group_averages_L1[name])

	# Take L2 average of each
	L2_pushed_arr = []
	for i in range(len(region_names)):
		group_arr = []
		for j in range(len(region_map[region_names[i]])):
			group_arr.append(np.array(sparse_matrix_L2[region_map[region_names[i]][j]].todense())[0])
		average = L2U.mean_of_vectors(group_arr)
		group_averages_L2[region_names[i]] = average
		L2_pushed_arr.append(average)

	# Store L2 averages
	print("\nL2 Group Averages:")
	CSV.write(output_file, ["L2 Group Averages:"])
	for name in region_names:
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {group_averages_L2[name]}")
		CSV.write(output_file, group_averages_L2[name])

	# Push L1 down and store
	print("\nL1 Inverse Push Up:")
	CSV.write(output_file, ["L1 Pushed Up:"])
	L1_neg_arr = []
	for name in region_names:
		neg_count = 0
		median_inverse = L1U.inverse_push_up(group_averages_L1[name], T1, l1, nodes_in_order)
		for i in range(len(median_inverse)):
			if median_inverse[i] < negatives_filtering_threshold:
				neg_count += 1
		L1_neg_arr.append(neg_count)
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {median_inverse}")
		CSV.write(output_file, median_inverse)

	# Push L2 down and store
	print("\nL2 Inverse Push Up:")
	CSV.write(output_file, ["L2 Pushed Up:"])
	L2_neg_arr = []
	for name in region_names:
		neg_count = 0
		mean_inverse = L2U.inverse_push_up(group_averages_L2[name], T1, l1, nodes_in_order)
		for i in range(len(mean_inverse)):
			if mean_inverse[i] < negatives_filtering_threshold:
				neg_count += 1
		L2_neg_arr.append(neg_count)
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {mean_inverse}")
		CSV.write(output_file, mean_inverse)

	# Write negative counts
	CSV.write(output_file, ["L1 and L2 Negatives by Group:"])
	CSV.write(output_file, L1_neg_arr)
	CSV.write(output_file, L2_neg_arr)

	L1_distance_matrix = compute_pairwise_pushed(L1_pushed_arr)
	L2_distance_matrix = compute_pairwise_pushed(L2_pushed_arr)

	print("L1 Distance Matrix:")
	CSV.write(output_file, ["L1 Distance Matrix:"])
	for i in range(len(L1_pushed_arr)):
		print(L1_distance_matrix[i])
		CSV.write(output_file, L1_distance_matrix[i])

	print("L2 Distance Matrix:")
	CSV.write(output_file, ["L2 Distance Matrix:"])
	for i in range(len(L2_pushed_arr)):
		print(L2_distance_matrix[i])
		CSV.write(output_file, L2_distance_matrix[i])		

# Argument parsing
if __name__ == "__main__":
	args = sys.argv
	if len(args) != 8:
		raise Exception("Invalid number of parameters.")
	else:
		L1_file = args[1]
		L2_file = args[2]
		biom_file = args[3]
		tree_file = args[4]
		metadata_file = args[5]
		tax_file = args[6]
		output_file = args[7]
		print(L1_file, L2_file, biom_file, tree_file, metadata_file, tax_file, output_file)
		if not path.exists(L1_file) or not path.exists(L2_file) or not path.exists(biom_file) or not path.exists(tree_file) or not path.exists(metadata_file) or not path.exists(tax_file):
			raise Exception("Error: Invalid file path(s).")
		compute_averages(L1_file, L2_file, biom_file, tree_file, metadata_file, tax_file, output_file)
