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

#nodes_samples = BW.extract_biom('../data/47422_otu_table.biom')
#T1, l1, nodes_in_order = L2U.parse_tree_file('../data/trees/gg_13_5_otus_99_annotated.tree')
#(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

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

def Group_Pairwise(debug):
	metadata = meta.extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt')
	sample_groups = []
	groups_temp = list(metadata.values())
	groups = []
	for i in range(len(groups_temp)):
		if groups_temp[i]['body_site'] not in groups:
			groups.append(groups_temp[i]['body_site'])
	print(groups)
	#for i in range(len(groups)):

	print(metadata, PCoA_Samples)

if __name__ == "__main__":

	args = sys.argv
	if len(args) == 3:
		if args[1] == "-v" or args[2] == "-v":
			debug = 1
		else:
			debug = 0
		if args[2] == "0" or args[1] == "0":
			toggle = 0
		else:
			toggle = 1
	elif len(args) == 2:
		if args[1] == "-v":
			debug = 1
		else:
			debug = 0
		if args[1] == "0":
			toggle = 0
		else:
			toggle = 1
	else:
		debug = 0
		toggle = 1

	if toggle:
		Group_Pairwise(debug)
	else:
		Total_Pairwise(debug)

	