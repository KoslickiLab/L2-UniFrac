import sys
sys.path.append('../src/')
from L2UniFrac import push_up
from data import extract_biom, extract_samples, write, parse_tree_file, parse_envs
import multiprocessing as mp

T1 = {}
l1 = {}
nodes_in_order = []
nodes_weighted = {}
samples = []

def L2_pushup_worker(sample_num):
	L2_Pushed = push_up(nodes_weighted[samples[sample_num]], T1, l1, nodes_in_order)
	return L2_Pushed

def generate_preprocessed(biom_file, tree_file, output_file=None, max_cores=int(mp.cpu_count()/4)):
	global T1
	global l1
	global nodes_in_order
	global nodes_weighted
	global samples

	if max_cores > mp.cpu_count() or max_cores <= 1:
		cores = mp.cpu_count()-1
	else:
		cores = max_cores

	nodes_samples = extract_biom(biom_file)
	T1, l1, nodes_in_order = parse_tree_file(tree_file)
	(nodes_weighted, samples_temp) = parse_envs(nodes_samples, nodes_in_order)

	samples = extract_samples(biom_file)

	# Multi Core Method
	L2_preprocessed = []

	values = range(len(samples))

	dim1 = len(samples)
	dim2 = len(L2_pushup_worker(0))

	if output_file is not None and unifrac_code == 0 or unifrac_code == 1:
		write(output_file, [dim1, dim2])
	L2_preprocessed.append([dim1, dim2])

	with mp.Pool(processes=cores) as pool:
		result = pool.map(L2_pushup_worker, values)

	for i in range(len(result)):
		for j in range(len(result[i])):
			if result[i][j] != 0:
				if output_file is not None:
					write(output_file, [i, j, result[i][j]])
				L2_preprocessed.append([i, j, result[i][j]])

	return L2_preprocessed

def generate_preprocessed_from_dict(biom_file, tree_file, sample_dict, output_file=None, max_cores=int(mp.cpu_count()/4)):
	global T1
	global l1
	global nodes_in_order
	global nodes_weighted
	global samples

	if max_cores > mp.cpu_count() or max_cores <= 1:
		cores = mp.cpu_count()-1
	else:
		cores = max_cores

	T1, l1, nodes_in_order = parse_tree_file(tree_file)
	nodes_weighted = class_dict

	samples = list(nodes_weighted.keys())

	# Multi Core Method
	L2_preprocessed = []

	values = range(len(samples))

	dim1 = len(samples)
	dim2 = len(L2_pushup_worker(0))

	if output_file is not None and unifrac_code == 0 or unifrac_code == 1:
		write(output_file, [dim1, dim2])
	L2_preprocessed.append([dim1, dim2])

	with mp.Pool(processes=cores) as pool:
		result = pool.map(L2_pushup_worker, values)

	for i in range(len(result)):
		for j in range(len(result[i])):
			if result[i][j] != 0:
				if output_file is not None:
					write(output_file, [i, j, result[i][j]])
				L2_preprocessed.append([i, j, result[i][j]])

	return L2_preprocessed

if __name__ == '__main__':
	generate_preprocessed('../data/47422_otu_table.biom', '../data/trees/gg_13_5_otus_99_annotated.tree', 'L2-Push-Out.csv')
	