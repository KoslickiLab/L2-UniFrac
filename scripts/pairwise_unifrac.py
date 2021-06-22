import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L2Unifrac as L2U
import BiomWrapper as BW
import write_to_csv as CSV
import metadata_wrapper as meta
import multiprocessing as mp

cores = mp.cpu_count()

nodes_samples = BW.extract_biom('../data/47422_otu_table.biom')
T1, l1, nodes_in_order = L2U.parse_tree_file('../data/trees/gg_13_5_otus_99_annotated.tree')
(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')

def unifrac_work_wrapper(args):
	return unifrac_worker(*args)

def unifrac_worker(samp1num, samp2num):
	L2UniFrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_in_order, nodes_weighted[PCoA_Samples[samp1num]], nodes_weighted[PCoA_Samples[samp2num]])
	formatted_L2 = "{:.16f}".format(L2UniFrac)
	return L2UniFrac, f"\tInner loop: {str(samp2num).zfill(4)} | L2-UniFrac: {formatted_L2} | Sample 1: {PCoA_Samples[samp1num]} | Sample 2: {PCoA_Samples[samp2num]}"

def Total_Pairwise(debug):
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

	for i in range(len(PCoA_Samples)):
		dist_list = []
		for j in range(len(PCoA_Samples)):
			dist_list.append(result[i*len(PCoA_Samples)+j][0])
			if debug == 1:
				print(result[i*len(PCoA_Samples)+j][1])

		CSV.write('L2-UniFrac-Out.csv', dist_list)

def Group_Pairwise(debug, group_num):
	group_num -= 1
	metadata = meta.extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt')
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

	for i in range(len(sample_sites[group_num])):
		dist_list = []
		for j in range(len(sample_sites[group_num])):
			dist_list.append(result[i*len(sample_sites[group_num])+j][0])
			if debug == 1:
				print(result[i*len(sample_sites[group_num])+j][1])

		CSV.write('L2-UniFrac-Out-' + groups[group_num] + '.csv', dist_list)

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
		Group_Pairwise(debug, group_num)
	else:
		pass
		Total_Pairwise(debug)

	