"""
This script aims to simply compute the a taxonomic profile of two pre-selected averages
in terms of their differential abundances in order to verify the differential abundance
computations.
"""

import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
sys.path.append('../scripts')
import L2Unifrac as L2U
import averages as avg
import TaxWrapper as tax
import numpy as np

# Compute diffab
(Tint, lint, nodes_in_order) = L2U.parse_tree_file('../data/trees/gg_13_5_otus_99_annotated.tree')
L2_region_names, L2_tax_arr, L2_group_averages, L2_inverse_pushed, L2_neg_arr, L2_distance_matrix, L2_node_type_group_abundances = avg.compute_L2_averages('../scripts/L2-Push-Out.csv', '../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt', '../data/taxonomies/gg_13_8_99.gg.tax', '../scripts/Group-Averages-2.csv', True)
L2_UniFrac, DifferentialAbundance = L2U.L2Unifrac_weighted(Tint, lint, nodes_in_order, L2_inverse_pushed[L2_region_names[0]], L2_inverse_pushed[L2_region_names[1]])

# Separate dictionaries
group_1_diffab = {key[0]:0 for key, value in DifferentialAbundance.items()}
group_2_diffab = {key[0]:0 for key, value in DifferentialAbundance.items()}
for key, value in DifferentialAbundance.items():
	if DifferentialAbundance[key] > 0:
		group_1_diffab[key[0]] = value
	elif DifferentialAbundance[key] < 0:
		group_2_diffab[key[0]] = -value

# Normalize to fractional amount
group_1_sum = sum(list(group_1_diffab.values()))
group_2_sum = sum(list(group_2_diffab.values()))
for key, value in group_1_diffab.items():
	group_1_diffab[key] = value/group_1_sum
for key, value in group_2_diffab.items():
	group_2_diffab[key] = value/group_2_sum

# Format taxonomies
new_L2_tax_arr = []
for i in range(len(L2_tax_arr)):
	new_L2_tax_arr.append('|'.join(L2_tax_arr[i].split(';')[:-1]))

taxonomies = tax.extract_tax('../data/taxonomies/gg_13_8_99.gg.tax')
tax_id_ref = {}
for key, value in taxonomies.items():
	tax_id_ref[value] = key

profile_list_1 = []
profile_list_1.append('@SampleID:sample_0')
profile_list_1.append('@Version:0.9')
profile_list_1.append('@Ranks: superkingdom|phylum|class|order|family|genus|species')
profile_list_1.append('')
profile_list_1.append('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE')
for key, value in group_1_diffab.items():
	split_tax = new_L2_tax_arr[key].split('|')
	outermost = 'Root'
	if split_tax[-1][0] == 'k':
		outermost = 'superkingdom'
	elif split_tax[-1][0] == 'p':
		outermost = 'phylum'
	elif split_tax[-1][0] == 'c':
		outermost = 'class'
	elif split_tax[-1][0] == 'o':
		outermost = 'order'
	elif split_tax[-1][0] == 'f':
		outermost = 'family'
	elif split_tax[-1][0] == 'g':
		outermost = 'genus'
	elif split_tax[-1][0] == 's':
		outermost = 'species'
	tax_path = ''
	tax_id = 0
	last_tax_id = 0
	for j in range(len(split_tax)):
		tmp_tax = ';'.join(split_tax[:j+1])+';'
		if len(tax_path) != 0:
			tax_path += '|'
		if tmp_tax in tax_id_ref:
			tax_path += str(tax_id_ref[tmp_tax])
			last_tax_id = str(tax_id_ref[tmp_tax])
		else:
			tax_path += last_tax_id

	profile_list_1.append('{0}\t{1}\t{2}\t{3}\t{4}'.format(last_tax_id, outermost, tax_path, new_L2_tax_arr[key], value))

with open('sample_0.profile', 'w') as f:
    for line in profile_list_1:
        f.write("{0}\n".format(line))

profile_list_2 = []
profile_list_2.append('@SampleID:sample_1')
profile_list_2.append('@Version:0.9')
profile_list_2.append('@Ranks: superkingdom|phylum|class|order|family|genus|species')
profile_list_2.append('')
profile_list_2.append('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE')
for key, value in group_2_diffab.items():
	split_tax = new_L2_tax_arr[key].split('|')
	outermost = 'Root'
	if split_tax[-1][0] == 'k':
		outermost = 'superkingdom'
	elif split_tax[-1][0] == 'p':
		outermost = 'phylum'
	elif split_tax[-1][0] == 'c':
		outermost = 'class'
	elif split_tax[-1][0] == 'o':
		outermost = 'order'
	elif split_tax[-1][0] == 'f':
		outermost = 'family'
	elif split_tax[-1][0] == 'g':
		outermost = 'genus'
	elif split_tax[-1][0] == 's':
		outermost = 'species'
	tax_path = ''
	tax_id = 0
	last_tax_id = 0
	for j in range(len(split_tax)):
		tmp_tax = ';'.join(split_tax[:j+1])+';'
		if len(tax_path) != 0:
			tax_path += '|'
		if tmp_tax in tax_id_ref:
			tax_path += str(tax_id_ref[tmp_tax])
			last_tax_id = str(tax_id_ref[tmp_tax])
		else:
			tax_path += last_tax_id

	profile_list_2.append('{0}\t{1}\t{2}\t{3}\t{4}'.format(last_tax_id, outermost, tax_path, new_L2_tax_arr[key], value))

with open('sample_1.profile', 'w') as f:
    for line in profile_list_2:
        f.write("{0}\n".format(line))