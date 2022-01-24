import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
from os import path
import os
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import CSVWrapper as CSV
import MetadataWrapper as meta
import TaxWrapper as tax
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix

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
negatives_filtering_threshold = -10e-12

# Check all descendents of temp node for most shared taxonomy
def most_shared_taxonomy(node, inverse_T1, nodes_in_order, taxonomies):
	if node not in inverse_T1:
		return taxonomies[int(nodes_in_order[node])]
	else:
		for i in range(len(inverse_T1[node])):
			pass

# Computes the pairwise average
def compute_pairwise_pushed_L1(pushed_arr):
	dist_matrix = []
	for i in range(len(pushed_arr)):
		dist_arr = []
		for j in range(len(pushed_arr)):
			unifrac_distance = np.sum(np.abs(pushed_arr[i] - pushed_arr[j]))
			dist_arr.append(unifrac_distance)
		dist_matrix.append(dist_arr)
	return dist_matrix

# Computes the pairwise average
def compute_pairwise_pushed_L2(pushed_arr):
	dist_matrix = []
	for i in range(len(pushed_arr)):
		dist_arr = []
		for j in range(len(pushed_arr)):
			unifrac_distance = np.linalg.norm(pushed_arr[i] - pushed_arr[j])
			dist_arr.append(unifrac_distance)
		dist_matrix.append(dist_arr)
	return dist_matrix

# Helper function for averages. Outputs a CSV containing info from each step and negative counts at end.
def compute_L1_averages(L1_file, biom_file, tree_file, metadata_file, tax_file, output_file=None, most_shared=False):
	
	if output_file is not None and path.exists(output_file):
		os.remove(output_file)

	# Note: these are the same for L1/L2, so they will be computed only once. (USE T1 FOR ANCESTORS FOR TEMP NODES)
	#nodes_samples = BW.extract_biom(biom_file)
	T1, l1, nodes_in_order = L1U.parse_tree_file(tree_file)

	# Subsample Biom file (2 samples). Then trace the mass to find where it no longer sums to 1.
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
	if not isinstance(L1_file, list):
		sparse_matrix_L1 = CSV.read_sparse(L1_file)
	else:
		sparse_matrix_L1 = L1_file

	group_averages_L1 = {}

	# Store region names for later
	if output_file is not None:
		CSV.write(output_file, region_names)

	# Write taxas for cell
	taxonomies = tax.extract_tax(tax_file)
	tax_arr = []

	if most_shared:
		leaf_nodes = []
		for i in range(len(nodes_in_order)):
			if nodes_in_order[i][0] != 't':
				leaf_nodes.append(i)
		descendent_dict = {i:[] for i in range(len(nodes_in_order))}
		for i in range(len(leaf_nodes)):
			temp_node = leaf_nodes[i]
			while True:
				if temp_node in T1:
					if temp_node != leaf_nodes[i]:
						descendent_dict[temp_node].append(leaf_nodes[i])
					temp_node = T1[temp_node]
				else:
					if temp_node != leaf_nodes[i]:
						descendent_dict[temp_node].append(leaf_nodes[i])
					break
		for i in range(len(nodes_in_order)):
			descendent_taxonomy = []
			for j in range(len(descendent_dict[i])):
				descendent_taxonomy.append(taxonomies[int(nodes_in_order[descendent_dict[i][j]])])
			if len(descendent_taxonomy) == 0:
				tax_arr.append(taxonomies[int(nodes_in_order[i])])
			else:
				k_dict = {}
				p_dict = {}
				c_dict = {}
				o_dict = {}
				f_dict = {}
				g_dict = {}
				s_dict = {}
				for j in range(len(descendent_taxonomy)):
					tmp_taxa = descendent_taxonomy[j].split(';')[:-1]
					for k in range(len(tmp_taxa)):
						if tmp_taxa[k][0] == 'k' and tmp_taxa[k] not in k_dict:
							k_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'k':
							k_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'p' and tmp_taxa[k] not in p_dict:
							p_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'p':
							p_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'c' and tmp_taxa[k] not in c_dict:
							c_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'c':
							c_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'o' and tmp_taxa[k] not in o_dict:
							o_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'o':
							o_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'f' and tmp_taxa[k] not in f_dict:
							f_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'f':
							f_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'g' and tmp_taxa[k] not in g_dict:
							g_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'g':
							g_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 's' and tmp_taxa[k] not in s_dict:
							s_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 's':
							s_dict[tmp_taxa[k]] += 1
				shared_taxonomy = ''
				for key, value in k_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+key
						break
				for key, value in p_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in c_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in o_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in f_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in g_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in s_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				shared_taxonomy = shared_taxonomy+';'
				if shared_taxonomy == ';':
					shared_taxonomy = 'Root' # Root node will include taxonomy from Archaea and Bacteria, thus sharing nothing.
				tax_arr.append(shared_taxonomy)
	else:
		for i in range(len(nodes_in_order)):
			if nodes_in_order[i][0] != 't':
				tax_arr.append(taxonomies[int(nodes_in_order[i])])
			else:
				loop = True
				if i in T1:
					temp_node = T1[i]
				else:
					tax_arr.append('internal')
					loop = False
				while loop:
					if nodes_in_order[temp_node][0] != 't':
						tax_arr.append(taxonomies[int(nodes_in_order[temp_node])])
						break
					else:
						if temp_node in T1:
							temp_node = T1[temp_node]
						else:
							tax_arr.append('internal')
							break
	
	if output_file is not None:
		CSV.write(output_file, tax_arr)

	# Take L1 average of each
	L1_pushed_arr = []
	for i in range(len(region_names)):
		group_arr = []
		if not isinstance(L1_file, list):
			for j in range(len(region_map[region_names[i]])):
				group_arr.append(np.array(sparse_matrix_L1[region_map[region_names[i]][j]].todense())[0])
		else:
			for j in range(len(region_map[region_names[i]])):
				group_arr.append(sparse_matrix_L1[region_map[region_names[i]][j]])
		average = L1U.median_of_vectors(group_arr)
		group_averages_L1[region_names[i]] = average
		L1_pushed_arr.append(average)

	# Store L1 averages
	print("L1 Group Averages:")
	if output_file is not None:
		CSV.write(output_file, ["L1 Group Averages:"])
	for name in region_names:
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {group_averages_L1[name]}")
		if output_file is not None:
			CSV.write(output_file, group_averages_L1[name])

	# Push L1 down and store
	print("\nL1 Inverse Push Up:")
	if output_file is not None:
		CSV.write(output_file, ["L1 Inverse Pushed Up:"])
	L1_neg_arr = []
	L1_inverse_pushed = {}
	for name in region_names:
		neg_count = 0
		median_inverse = L1U.inverse_push_up(group_averages_L1[name], T1, l1, nodes_in_order)
		L1_inverse_pushed[name] = median_inverse
		for i in range(len(median_inverse)):
			if median_inverse[i] < negatives_filtering_threshold:
				neg_count += 1
		L1_neg_arr.append(neg_count)
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {median_inverse}")
		if output_file is not None:
			CSV.write(output_file, median_inverse)

	# Write negative counts
	print("L1 Negatives by Group:")
	print(L1_neg_arr)
	if output_file is not None:
		CSV.write(output_file, ["L1 Negatives by Group:"])
		CSV.write(output_file, L1_neg_arr)

	L1_distance_matrix = compute_pairwise_pushed_L1(L1_pushed_arr)

	print("L1 Distance Matrix:")
	if output_file is not None:
		CSV.write(output_file, ["L1 Distance Matrix:"])
	for i in range(len(L1_pushed_arr)):
		print(L1_distance_matrix[i])
		if output_file is not None:
			CSV.write(output_file, L1_distance_matrix[i])	

	print("L1 Abundances by Node Type:")
	if output_file is not None:
		CSV.write(output_file, ["L1 Abundances by Node Type:"])
	L1_node_type_group_abundances = []
	for name in region_names:
		region_abundance_vector = L1_inverse_pushed[name]
		k = p = c = o = f = g = s = temp = 0
		for i in range(len(region_abundance_vector)):
			node_tax = tax_arr[i].split(';')
			if len(node_tax) > 1:
				if node_tax[-2][0] == 'k':
					k += region_abundance_vector[i]
				elif node_tax[-2][0] == 'p':
					p += region_abundance_vector[i]
				elif node_tax[-2][0] == 'c':
					c += region_abundance_vector[i]
				elif node_tax[-2][0] == 'o':
					o += region_abundance_vector[i]
				elif node_tax[-2][0] == 'f':
					f += region_abundance_vector[i]
				elif node_tax[-2][0] == 'g':
					g += region_abundance_vector[i]
				elif node_tax[-2][0] == 's':
					s += region_abundance_vector[i]
				else:
					print("Error")
			else:
				temp += region_abundance_vector[i]
		print([k, p, c, o, f, g, s, temp])
		if output_file is not None:
			CSV.write(output_file, [k, p, c, o, f, g, s, temp])
		L1_node_type_group_abundances.append([k, p, c, o, f, g, s, temp])

	return region_names, tax_arr, group_averages_L1, L1_inverse_pushed, L1_neg_arr, L1_distance_matrix, L1_node_type_group_abundances

# Helper function for averages. Outputs a CSV containing info from each step and negative counts at end.
def compute_L2_averages(L2_file, biom_file, tree_file, metadata_file, tax_file, output_file=None, most_shared=False):

	if output_file is not None and path.exists(output_file):
		os.remove(output_file)

	# Note: these are the same for L1/L2, so they will be computed only once. (USE T1 FOR ANCESTORS FOR TEMP NODES)
	#nodes_samples = BW.extract_biom(biom_file)
	T1, l1, nodes_in_order = L2U.parse_tree_file(tree_file)
	# Subsample Biom file (2 samples). Then trace the mass to find where it no longer sums to 1.
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
	# Read sparse matrix
	if not isinstance(L2_file, list):
		sparse_matrix_L2 = CSV.read_sparse(L2_file)
	else:
		sparse_matrix_L2 = L2_file

	group_averages_L2 = {}
	# Store region names for later
	if output_file is not None:
		CSV.write(output_file, region_names)
	# Write taxas for cell
	taxonomies = tax.extract_tax(tax_file)
	tax_arr = []

	if most_shared:
		leaf_nodes = []
		for i in range(len(nodes_in_order)):
			if nodes_in_order[i][0] != 't':
				leaf_nodes.append(i)
		descendent_dict = {i:[] for i in range(len(nodes_in_order))}
		for i in range(len(leaf_nodes)):
			temp_node = leaf_nodes[i]
			while True:
				if temp_node in T1:
					if temp_node != leaf_nodes[i]:
						descendent_dict[temp_node].append(leaf_nodes[i])
					temp_node = T1[temp_node]
				else:
					if temp_node != leaf_nodes[i]:
						descendent_dict[temp_node].append(leaf_nodes[i])
					break
		for i in range(len(nodes_in_order)):
			descendent_taxonomy = []
			for j in range(len(descendent_dict[i])):
				descendent_taxonomy.append(taxonomies[int(nodes_in_order[descendent_dict[i][j]])])
			if len(descendent_taxonomy) == 0:
				tax_arr.append(taxonomies[int(nodes_in_order[i])])
			else:
				k_dict = {}
				p_dict = {}
				c_dict = {}
				o_dict = {}
				f_dict = {}
				g_dict = {}
				s_dict = {}
				for j in range(len(descendent_taxonomy)):
					tmp_taxa = descendent_taxonomy[j].split(';')[:-1]
					for k in range(len(tmp_taxa)):
						if tmp_taxa[k][0] == 'k' and tmp_taxa[k] not in k_dict:
							k_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'k':
							k_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'p' and tmp_taxa[k] not in p_dict:
							p_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'p':
							p_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'c' and tmp_taxa[k] not in c_dict:
							c_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'c':
							c_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'o' and tmp_taxa[k] not in o_dict:
							o_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'o':
							o_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'f' and tmp_taxa[k] not in f_dict:
							f_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'f':
							f_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 'g' and tmp_taxa[k] not in g_dict:
							g_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 'g':
							g_dict[tmp_taxa[k]] += 1
						if tmp_taxa[k][0] == 's' and tmp_taxa[k] not in s_dict:
							s_dict[tmp_taxa[k]] = 1
						elif tmp_taxa[k][0] == 's':
							s_dict[tmp_taxa[k]] += 1
				shared_taxonomy = ''
				for key, value in k_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+key
						break
				for key, value in p_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in c_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in o_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in f_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in g_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				for key, value in s_dict.items():
					if value == len(descendent_taxonomy):
						shared_taxonomy = shared_taxonomy+';'+key
						break
				shared_taxonomy = shared_taxonomy+';'
				if shared_taxonomy == ';':
					shared_taxonomy = 'Root' # Root node will include taxonomy from Archaea and Bacteria, thus sharing nothing.
				tax_arr.append(shared_taxonomy)
	else:
		for i in range(len(nodes_in_order)):
			if nodes_in_order[i][0] != 't':
				tax_arr.append(taxonomies[int(nodes_in_order[i])])
			else:
				loop = True
				if i in T1:
					temp_node = T1[i]
				else:
					tax_arr.append('internal')
					loop = False
				while loop:
					if nodes_in_order[temp_node][0] != 't':
						tax_arr.append(taxonomies[int(nodes_in_order[temp_node])])
						break
					else:
						if temp_node in T1:
							temp_node = T1[temp_node]
						else:
							tax_arr.append('internal')
							break
	return

	if output_file is not None:
		CSV.write(output_file, tax_arr)

	# Take L2 average of each
	L2_pushed_arr = []
	for i in range(len(region_names)):
		group_arr = []
		if not isinstance(L2_file, list):
			for j in range(len(region_map[region_names[i]])):
				group_arr.append(np.array(sparse_matrix_L2[region_map[region_names[i]][j]].todense())[0])
		else:
			for j in range(len(region_map[region_names[i]])):
				group_arr.append(np.array(sparse_matrix_L2[region_map[region_names[i]][j]])[0])
		average = L2U.mean_of_vectors(group_arr)
		group_averages_L2[region_names[i]] = average
		L2_pushed_arr.append(average)

	# Store L2 averages
	print("\nL2 Group Averages:")
	if output_file is not None:
		CSV.write(output_file, ["L2 Group Averages:"])
	for name in region_names:
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {group_averages_L2[name]}")
		if output_file is not None:
			CSV.write(output_file, group_averages_L2[name])

	# Push L2 down and store
	print("\nL2 Inverse Push Up:")
	if output_file is not None:
		CSV.write(output_file, ["L2 Inverse Pushed Up:"])
	L2_neg_arr = []
	L2_inverse_pushed = {}
	for name in region_names:
		neg_count = 0
		mean_inverse = L2U.inverse_push_up(group_averages_L2[name], T1, l1, nodes_in_order)
		L2_inverse_pushed[name] = mean_inverse
		for i in range(len(mean_inverse)):
			if mean_inverse[i] < negatives_filtering_threshold:
				neg_count += 1
		L2_neg_arr.append(neg_count)
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {mean_inverse}")
		if output_file is not None:
			CSV.write(output_file, mean_inverse)

	# Write negative counts
	print("L2 Negatives by Group:")
	print(L2_neg_arr)
	if output_file is not None:
		CSV.write(output_file, ["L2 Negatives by Group:"])
		CSV.write(output_file, L2_neg_arr)

	L2_distance_matrix = compute_pairwise_pushed_L2(L2_pushed_arr)

	print("L2 Distance Matrix:")
	if output_file is not None:
		CSV.write(output_file, ["L2 Distance Matrix:"])
	for i in range(len(L2_pushed_arr)):
		print(L2_distance_matrix[i])
		if output_file is not None:
			CSV.write(output_file, L2_distance_matrix[i])		

	print("L2 Abundances by Node Type:")
	if output_file is not None:
		CSV.write(output_file, ["L2 Abundances by Node Type:"])
	L2_node_type_group_abundances = []
	for name in region_names:
		region_abundance_vector = L2_inverse_pushed[name]
		k = p = c = o = f = g = s = temp = 0
		for i in range(len(region_abundance_vector)):
			node_tax = tax_arr[i].split(';')
			if len(node_tax) > 1:
				if node_tax[-2][0] == 'k':
					k += region_abundance_vector[i]
				elif node_tax[-2][0] == 'p':
					p += region_abundance_vector[i]
				elif node_tax[-2][0] == 'c':
					c += region_abundance_vector[i]
				elif node_tax[-2][0] == 'o':
					o += region_abundance_vector[i]
				elif node_tax[-2][0] == 'f':
					f += region_abundance_vector[i]
				elif node_tax[-2][0] == 'g':
					g += region_abundance_vector[i]
				elif node_tax[-2][0] == 's':
					s += region_abundance_vector[i]
				else:
					print("Error")
			else:
				temp += region_abundance_vector[i]
		print([k, p, c, o, f, g, s, temp])
		if output_file is not None:
			CSV.write(output_file, [k, p, c, o, f, g, s, temp])
		L2_node_type_group_abundances.append([k, p, c, o, f, g, s, temp])

	return region_names, tax_arr, group_averages_L2, L2_inverse_pushed, L2_neg_arr, L2_distance_matrix, L2_node_type_group_abundances

# Helper function for averages. Outputs a CSV containing info from each step and negative counts at end.
def compute_L1_L2_averages(L1_file, L2_file, biom_file, tree_file, metadata_file, tax_file, output_file=None):

	if output_file is not None and path.exists(output_file):
		os.remove(output_file)

	# Note: these are the same for L1/L2, so they will be computed only once. (USE T1 FOR ANCESTORS FOR TEMP NODES)
	#nodes_samples = BW.extract_biom(biom_file)
	T1, l1, nodes_in_order = L2U.parse_tree_file(tree_file)

	# Subsample Biom file (2 samples). Then trace the mass to find where it no longer sums to 1.
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
	if not isinstance(L1_file, list):
		sparse_matrix_L1 = CSV.read_sparse(L1_file)
	else:
		sparse_matrix_L1 = L1_file
	if not isinstance(L2_file, list):
		sparse_matrix_L2 = CSV.read_sparse(L2_file)
	else:
		sparse_matrix_L2 = L2_file

	group_averages_L1 = {}
	group_averages_L2 = {}

	# Store region names for later
	if output_file is not None:
		CSV.write(output_file, region_names)

	# Write taxas for cell
	taxonomies = tax.extract_tax(tax_file)
	tax_arr = []
	for i in range(len(nodes_in_order)):
		if nodes_in_order[i][0] != 't':
			tax_arr.append(taxonomies[int(nodes_in_order[i])])
		else:
			loop = True
			if i in T1:
				temp_node = T1[i]
			else:
				tax_arr.append('internal')
				loop = False
			while loop:
				if nodes_in_order[temp_node][0] != 't':
					tax_arr.append(taxonomies[int(nodes_in_order[temp_node])])
					break
				else:
					if temp_node in T1:
						temp_node = T1[temp_node]
					else:
						tax_arr.append('internal')
						break
	
	if output_file is not None:
		CSV.write(output_file, tax_arr)

	# Take L1 average of each
	L1_pushed_arr = []
	for i in range(len(region_names)):
		group_arr = []
		if not isinstance(L1_file, list):
			for j in range(len(region_map[region_names[i]])):
				group_arr.append(np.array(sparse_matrix_L1[region_map[region_names[i]][j]].todense())[0])
		else:
			for j in range(len(region_map[region_names[i]])):
				group_arr.append(np.array(sparse_matrix_L1[region_map[region_names[i]][j]])[0])
		average = L1U.median_of_vectors(group_arr)
		group_averages_L1[region_names[i]] = average
		L1_pushed_arr.append(average)

	# Store L1 averages
	print("L1 Group Averages:")
	if output_file is not None:
		CSV.write(output_file, ["L1 Group Averages:"])
	for name in region_names:
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {group_averages_L1[name]}")
		if output_file is not None:
			CSV.write(output_file, group_averages_L1[name])

	# Take L2 average of each
	L2_pushed_arr = []
	for i in range(len(region_names)):
		group_arr = []
		if not isinstance(L2_file, list):
			for j in range(len(region_map[region_names[i]])):
				group_arr.append(np.array(sparse_matrix_L2[region_map[region_names[i]][j]].todense())[0])
		else:
			for j in range(len(region_map[region_names[i]])):
				group_arr.append(np.array(sparse_matrix_L2[region_map[region_names[i]][j]])[0])
		average = L2U.mean_of_vectors(group_arr)
		group_averages_L2[region_names[i]] = average
		L2_pushed_arr.append(average)

	# Store L2 averages
	print("\nL2 Group Averages:")
	if output_file is not None:
		CSV.write(output_file, ["L2 Group Averages:"])
	for name in region_names:
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {group_averages_L2[name]}")
		if output_file is not None:
			CSV.write(output_file, group_averages_L2[name])

	# Push L1 down and store
	print("\nL1 Inverse Push Up:")
	if output_file is not None:
		CSV.write(output_file, ["L1 Inverse Pushed Up:"])
	L1_neg_arr = []
	L1_inverse_pushed = {}
	for name in region_names:
		neg_count = 0
		median_inverse = L1U.inverse_push_up(group_averages_L1[name], T1, l1, nodes_in_order)
		L1_inverse_pushed[name] = median_inverse
		for i in range(len(median_inverse)):
			if median_inverse[i] < negatives_filtering_threshold:
				neg_count += 1
		L1_neg_arr.append(neg_count)
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {median_inverse}")
		if output_file is not None:
			CSV.write(output_file, median_inverse)

	# Push L2 down and store
	print("\nL2 Inverse Push Up:")
	if output_file is not None:
		CSV.write(output_file, ["L2 Inverse Pushed Up:"])
	L2_neg_arr = []
	L2_inverse_pushed = {}
	for name in region_names:
		neg_count = 0
		mean_inverse = L2U.inverse_push_up(group_averages_L2[name], T1, l1, nodes_in_order)
		L2_inverse_pushed[name] = mean_inverse
		for i in range(len(mean_inverse)):
			if mean_inverse[i] < negatives_filtering_threshold:
				neg_count += 1
		L2_neg_arr.append(neg_count)
		padded_name = "{:<15}".format(name+":")
		print(f"{padded_name} {mean_inverse}")
		if output_file is not None:
			CSV.write(output_file, mean_inverse)

	# Write negative counts
	print("L1 and L2 Negatives by Group:")
	print(L1_neg_arr)
	print(L2_neg_arr)
	if output_file is not None:
		CSV.write(output_file, ["L1 and L2 Negatives by Group:"])
		CSV.write(output_file, L1_neg_arr)
		CSV.write(output_file, L2_neg_arr)

	L1_distance_matrix = compute_pairwise_pushed_L1(L1_pushed_arr)
	L2_distance_matrix = compute_pairwise_pushed_L2(L2_pushed_arr)

	print("L1 Distance Matrix:")
	if output_file is not None:
		CSV.write(output_file, ["L1 Distance Matrix:"])
	for i in range(len(L1_pushed_arr)):
		print(L1_distance_matrix[i])
		if output_file is not None:
			CSV.write(output_file, L1_distance_matrix[i])

	print("L2 Distance Matrix:")
	if output_file is not None:
		CSV.write(output_file, ["L2 Distance Matrix:"])
	for i in range(len(L2_pushed_arr)):
		print(L2_distance_matrix[i])
		if output_file is not None:
			CSV.write(output_file, L2_distance_matrix[i])		

	print("L1 Abundances by Node Type:")
	if output_file is not None:
		CSV.write(output_file, ["L1 Abundances by Node Type:"])
	L1_node_type_group_abundances = []
	for name in region_names:
		region_abundance_vector = L1_inverse_pushed[name]
		k = p = c = o = f = g = s = temp = 0
		for i in range(len(region_abundance_vector)):
			node_tax = tax_arr[i].split(';')
			if len(node_tax) > 1:
				if node_tax[-2][0] == 'k':
					k += region_abundance_vector[i]
				elif node_tax[-2][0] == 'p':
					p += region_abundance_vector[i]
				elif node_tax[-2][0] == 'c':
					c += region_abundance_vector[i]
				elif node_tax[-2][0] == 'o':
					o += region_abundance_vector[i]
				elif node_tax[-2][0] == 'f':
					f += region_abundance_vector[i]
				elif node_tax[-2][0] == 'g':
					g += region_abundance_vector[i]
				elif node_tax[-2][0] == 's':
					s += region_abundance_vector[i]
				else:
					print("Error")
			else:
				temp += region_abundance_vector[i]
		print([k, p, c, o, f, g, s, temp])
		if output_file is not None:
			CSV.write(output_file, [k, p, c, o, f, g, s, temp])
		L1_node_type_group_abundances.append([k, p, c, o, f, g, s, temp])

	print("L2 Abundances by Node Type:")
	if output_file is not None:
		CSV.write(output_file, ["L2 Abundances by Node Type:"])
	L2_node_type_group_abundances = []
	for name in region_names:
		region_abundance_vector = L2_inverse_pushed[name]
		k = p = c = o = f = g = s = temp = 0
		for i in range(len(region_abundance_vector)):
			node_tax = tax_arr[i].split(';')
			if len(node_tax) > 1:
				if node_tax[-2][0] == 'k':
					k += region_abundance_vector[i]
				elif node_tax[-2][0] == 'p':
					p += region_abundance_vector[i]
				elif node_tax[-2][0] == 'c':
					c += region_abundance_vector[i]
				elif node_tax[-2][0] == 'o':
					o += region_abundance_vector[i]
				elif node_tax[-2][0] == 'f':
					f += region_abundance_vector[i]
				elif node_tax[-2][0] == 'g':
					g += region_abundance_vector[i]
				elif node_tax[-2][0] == 's':
					s += region_abundance_vector[i]
				else:
					print("Error")
			else:
				temp += region_abundance_vector[i]
		print([k, p, c, o, f, g, s, temp])
		if output_file is not None:
			CSV.write(output_file, [k, p, c, o, f, g, s, temp])
		L2_node_type_group_abundances.append([k, p, c, o, f, g, s, temp])

	return region_names, tax_arr, group_averages_L1, group_averages_L2, L1_inverse_pushed, L2_inverse_pushed, L1_neg_arr, L2_neg_arr, L1_distance_matrix, L2_distance_matrix, L1_node_type_group_abundances, L2_node_type_group_abundances

# Argument parsing
if __name__ == "__main__":
	#compute_L1_L2_averages('L1-Push-Out.csv', 'L2-Push-Out.csv', '../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt', '../data/taxonomies/gg_13_8_99.gg.tax', 'Group-Averages.csv')
	compute_L2_averages('L2-Push-Out.csv', '../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt', '../data/taxonomies/gg_13_8_99.gg.tax', 'Group-Averages-2.csv', True)
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
		compute_L1_L2_averages(L1_file, L2_file, biom_file, tree_file, metadata_file, tax_file, output_file)
