import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L2Unifrac as L2U
import BiomWrapper as BW
import CSVWrapper as CSV
import MetadataWrapper as meta
import multiprocessing as mp

T1 = {}
l1 = {}
nodes_in_order = []
nodes_weighted = {}

def unifrac_work_wrapper(args):
	return unifrac_worker(*args)

def unifrac_worker(samp1num, samp2num):
	L2UniFrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_in_order, nodes_weighted[PCoA_Samples[samp1num]], nodes_weighted[PCoA_Samples[samp2num]])
	formatted_L2 = "{:.16f}".format(L2UniFrac)
	return L2UniFrac, f"\tInner loop: {str(samp2num).zfill(4)} | L2-UniFrac: {formatted_L2} | Sample 1: {PCoA_Samples[samp1num]} | Sample 2: {PCoA_Samples[samp2num]}"

def Total_Pairwise(biom_file, tree_file, output_file=None, debug=0):
	global T1
	global l1
	global nodes_in_order
	global nodes_weighted

	cores = mp.cpu_count()

	nodes_samples = BW.extract_biom(biom_file)
	T1, l1, nodes_in_order = L2U.parse_tree_file(tree_file)
	(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

	PCoA_Samples = BW.extract_samples(biom_file)

	if debug == 1:
		print(f"Running Debugging Multiprocess on {cores-1} Cores...")

		# Testing subset of samples...
		PCoA_Samples = PCoA_Samples[:64]
		
		local_vars = list(locals().items())
		for var, obj in local_vars:
			print(f"{var.ljust(17)}: {sys.getsizeof(obj)}")

	# Multi Core Method
	row = [(i, j) for j in range(len(PCoA_Samples)) for i in range(len(PCoA_Samples))]

	with mp.Pool(processes=cores-1) as pool:
		result = pool.map(unifrac_work_wrapper, row)

	result_matrix = []
	for i in range(len(PCoA_Samples)):
		dist_list = []
		for j in range(len(PCoA_Samples)):
			dist_list.append(result[i*len(PCoA_Samples)+j][0])
			if debug == 1:
				print(result[i*len(PCoA_Samples)+j][1])
		result_matrix.append(dist_list)
		if output_file is not None:
			CSV.write(output_file, dist_list)
	return result_matrix

def Group_Pairwise(biom_file, tree_file, metadata_file, group_num, output_file=None, debug=0):
	global T1
	global l1
	global nodes_in_order
	global nodes_weighted

	cores = mp.cpu_count()

	nodes_samples = BW.extract_biom(biom_file)
	T1, l1, nodes_in_order = L2U.parse_tree_file(tree_file)
	(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

	PCoA_Samples = BW.extract_samples(biom_file)

	group_num -= 1
	metadata = meta.extract_metadata(metadata_file)
	sample_groups = []
	groups_temp = list(metadata.values())
	groups = []
	for i in range(len(groups_temp)):
		if groups_temp[i]['body_site'] not in groups:
			groups.append(groups_temp[i]['body_site'])
	print(groups)

	sample_sites = [[] for i in range(len(groups))]

	# Separate the groups
	for i in range(len(PCoA_Samples)):
		for j in range(len(groups)):
			if metadata[PCoA_Samples[i]]['body_site'] == groups[j]:
				sample_sites[j].append(PCoA_Samples[i])
	print(sample_sites)

	if debug == 1:
		print(f"Running Debugging Multiprocess on {cores-1} Cores...")

		# Testing subset of samples...
		sample_sites[group_num] = sample_sites[group_num][:64]
		
		local_vars = list(locals().items())
		for var, obj in local_vars:
			print(f"{var.ljust(17)}: {sys.getsizeof(obj)}")

	# Multi Core Method
	row = [(i, j) for j in range(len(sample_sites[group_num])) for i in range(len(sample_sites[group_num]))]

	with mp.Pool(processes=cores-1) as pool:
		result = pool.map(unifrac_work_wrapper, row)

	result_matrix = []
	for i in range(len(sample_sites[group_num])):
		dist_list = []
		for j in range(len(sample_sites[group_num])):
			dist_list.append(result[i*len(sample_sites[group_num])+j][0])
			if debug == 1:
				print(result[i*len(sample_sites[group_num])+j][1])
		result_matrix.append(dist_list)
		if output_file is not None:
			CSV.write(output_file[:-4] + "-"+ groups[group_num] + '.csv', dist_list)
	return result_matrix

if __name__ == "__main__":

	args = sys.argv
	for i in range(len(args)):
		try:
			args[i] = int(args[i])
		except:
			pass
	if len(args) > 3:
		raise Exception("Invalid number of parameters.")
	
	if "-h" in args:
		debug = 1
	else:
		debug = 0

	if len(args) > 1 and isinstance(args[1], int):
		group_num = args[1]
	elif len(args) > 2 and isinstance(args[2], int):
		group_num = args[2]
	else:
		group_num = 0

	if group_num:
		Group_Pairwise('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', '../data/metadata/P_1928_65684500_raw_meta.txt', group_num, 'L2-UniFrac-Out.csv', debug)
	else:
		Total_Pairwise('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', 'L2-UniFrac-Out.csv', debug)

	