import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import CSVWrapper as CSV
import multiprocessing as mp

def L1_pushup_worker(sample_num):
	L1_Pushed = L1U.push_up(nodes_weighted[PCoA_Samples[sample_num]], T1, l1, nodes_in_order)
	return L1_Pushed

def L2_pushup_worker(sample_num):
	L2_Pushed = L2U.push_up(nodes_weighted[PCoA_Samples[sample_num]], T1, l1, nodes_in_order)
	return L2_Pushed

def generate_preprocessed(biom_file, tree_file, unifrac_code, output_file_L1=None, output_file_L2=None, max_cores=int(mp.cpu_count()/4)):
	if max_cores > mp.cpu_count() or max_cores <= 1:
		cores = mp.cpu_count()-1
	else:
		cores = max_cores

	# Note: these are the same for L1/L2, so they will be computed only once.
	nodes_samples = BW.extract_biom(biom_file)
	T1, l1, nodes_in_order = L2U.parse_tree_file(tree_file)
	(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

	PCoA_Samples = BW.extract_samples(biom_file)

	# Multi Core Method
	L1_preprocessed = []
	L2_preprocessed = []

	values = range(len(PCoA_Samples))

	dim1 = len(PCoA_Samples)
	dim2 = len(L1_pushup_worker(0))

	if output_file_L1 is not None and unifrac_code == 1 or unifrac_code == 2:
		CSV.write(output_file_L1, [dim1, dim2])
	L1_preprocessed.append([dim1, dim2])
	if output_file_L2 is not None and unifrac_code == 0 or unifrac_code == 1:
		CSV.write(output_file_L2, [dim1, dim2])
	L2_preprocessed.append([dim1, dim2])

	if unifrac_code == 1 or unifrac_code == 2:

		with mp.Pool(processes=cores) as pool:
			result = pool.map(L1_pushup_worker, values)

		for i in range(len(result)):
			for j in range(len(result[i])):
				if result[i][j] != 0:
					if output_file_L1 is not None:
						CSV.write(output_file_L1, [i, j, result[i][j]])
					L1_preprocessed.append([i, j, result[i][j]])

	if unifrac_code == 0 or unifrac_code == 1:

		with mp.Pool(processes=cores) as pool:
			result = pool.map(L2_pushup_worker, values)

		for i in range(len(result)):
			for j in range(len(result[i])):
				if result[i][j] != 0:
					if output_file_L2 is not None:
						CSV.write(output_file_L2, [i, j, result[i][j]])
					L2_preprocessed.append([i, j, result[i][j]])

	if unifrac_code == 0:
		return [], L2_preprocessed
	if unifrac_code == 2:
		return L1_preprocessed, []
	if unifrac_code == 1:
		return L1_preprocessed, L2_preprocessed


def generate_preprocessed_combined(biom_file, tree_file, output_file_L1=None, output_file_L2=None, max_cores=int(mp.cpu_count()/4)):
	if max_cores > mp.cpu_count() or max_cores <= 1:
		cores = mp.cpu_count()-1
	else:
		cores = max_cores

	# Note: these are the same for L1/L2, so they will be computed only once.
	nodes_samples = BW.extract_biom(biom_file)
	T1, l1, nodes_in_order = L2U.parse_tree_file(tree_file)
	(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

	PCoA_Samples = BW.extract_samples(biom_file)

	# Multi Core Method
	L1_preprocessed = []
	L2_preprocessed = []

	values = range(len(PCoA_Samples))

	dim1 = len(PCoA_Samples)
	dim2 = len(L1_pushup_worker(0))

	if output_file_L1 is not None:
		CSV.write(output_file_L1, [dim1, dim2])
	L1_preprocessed.append([dim1, dim2])
	if output_file_L2 is not None:
		CSV.write(output_file_L2, [dim1, dim2])
	L2_preprocessed.append([dim1, dim2])

	with mp.Pool(processes=cores) as pool:
		L1_result = pool.map(L1_pushup_worker, values)

	for i in range(len(L1_result)):
		for j in range(len(L1_result[i])):
			if L1_result[i][j] != 0:
				if output_file_L1 is not None:
					CSV.write(output_file_L1, [i, j, L1_result[i][j]])
				L1_preprocessed.append([i, j, L1_result[i][j]])

	with mp.Pool(processes=int(cores/4-1)) as pool:
		L2_result = pool.map(L2_pushup_worker, values)

	for i in range(len(L2_result)):
		for j in range(len(L2_result[i])):
			if L2_result[i][j] != 0:
				if output_file_L2 is not None:
					CSV.write(output_file_L2, [i, j, L2_result[i][j]])
				L2_preprocessed.append([i, j, L2_result[i][j]])

	return L1_preprocessed, L2_preprocessed

if __name__ == '__main__':
	generate_preprocessed_combined('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', 'L1-Push-Out.csv', 'L2-Push-Out.csv')
	