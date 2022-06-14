import biom, csv, dendropy
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix

def extract_biom(address):

	# Load table into the biom Table format.
	Biom = biom.load_table(address)

	# Grab the sample IDs (i.e. '1928.SRS045191.SRX021470.SRR052699')
	sample_ids = Biom.ids()

	# Grab the node IDs for the tree (i.e. '858026')
	phylogenetic_tree_nodes = Biom.ids(axis='observation')

	# Formulate the dictionary required for L2 Unifrac by associating each sample weight to its node ID. 
	nodes_samples = {}
	nodes_samples_temp = {}
	for i in range(len(phylogenetic_tree_nodes)):
		row_vector = Biom._get_row(i).todense().tolist()[0] # Convert Table row to a standard python list.
		nodes_samples[phylogenetic_tree_nodes[i]] = {sample_ids[j]:row_vector[j] for j in range(len(sample_ids))}
		#sample_count = sum(nodes_samples_temp[phylogenetic_tree_nodes[i]].values())
		#nodes_samples[phylogenetic_tree_nodes[i]] = {k:v/sample_count for k,v in nodes_samples_temp[phylogenetic_tree_nodes[i]].items()} # NORMALIZE COLUMN, NOT ROW
	totals = [0 for i in range(len(nodes_samples[phylogenetic_tree_nodes[0]].keys()))]
	for i in range(len(phylogenetic_tree_nodes)):
		row = nodes_samples[phylogenetic_tree_nodes[i]]
		for index, (key, value) in enumerate(row.items()):
			totals[index] += value
	for i in range(len(phylogenetic_tree_nodes)):
		row = nodes_samples[phylogenetic_tree_nodes[i]]
		for index, (key, value) in enumerate(row.items()):
			nodes_samples[phylogenetic_tree_nodes[i]][key] /= totals[index]
	return nodes_samples

def extract_samples(address):

	# Load table into the biom Table format.
	Biom = biom.load_table(address)

	# Grab the sample IDs (i.e. '1928.SRS045191.SRX021470.SRR052699')
	sample_ids = Biom.ids()

	return sample_ids

def extract_nodes(address):

	# Load table into the biom Table format.
	Biom = biom.load_table(address)

	# Grab the node IDs for the tree (i.e. '858026')
	phylogenetic_tree_nodes = Biom.ids(axis='observation')

	return phylogenetic_tree_nodes

def extract_tax(extension):
	try:
		with open(str(extension)) as f:
			taxonomies = {}
			tax_list = f.read().splitlines() 
			for tax in range(len(tax_list)):
				tax_list[tax] = tax_list[tax].split('\t')
				taxonomies[int(tax_list[tax][0])] = tax_list[tax][1]
			return taxonomies
	except:
		raise Exception("Unknown exception occurred. Try again.")

# Function to pull metadata from specific metadata file sourced at http://mse.ac.cn/index.php/mse
def extract_metadata(extension):
	try:
		f = open(str(extension), 'r')
		metadata = f.readlines()

		for i in range(len(metadata)):
			metadata[i] = metadata[i].split("\t")

		del[metadata[0]]

		meta_dict = {}
		for i in range(len(metadata)):
			meta_dict[metadata[i][0]] = {'body_site': metadata[i][3][7:]}

		return meta_dict

	except FileNotFoundError:
		raise FileNotFoundError("Unknown file. Make sure you have the correct address and try again!")

	except:
		raise Exception("Unknown exception occurred. Try again.")

def extract_num_clusters(extension):
	try:
		f = open(str(extension), 'r')
		metadata = f.readlines()

		for i in range(len(metadata)):
			metadata[i] = metadata[i].split("\t")

		del[metadata[0]]

		meta_dict = {}
		for i in range(len(metadata)):
			meta_dict[metadata[i][0]] = {'body_site': metadata[i][3][7:]}

		values = list(meta_dict.values())
		
		unique_list = []
		for i in range(len(values)):
			if values[i]['body_site'] not in unique_list:
				unique_list.append(values[i]['body_site'])
		
		return len(unique_list)

	except FileNotFoundError:
		raise FileNotFoundError("Unknown file. Make sure you have the correct address and try again!")

	except:
		raise Exception("Unknown exception occurred. Try again.")

# Ensures that extracted metadata is suffient for the selected biom file
def test_metadata_completeness(metadata_path, biom_extension):
	meta_dict = extract_metadata(metadata_path)
	samples = extract_samples(biom_extension)

	for i in range(len(samples)):
		assert samples[i] in meta_dict

def write(name, dist_list):
	try:
		with open(name, 'a', newline='') as csvfile:
			file_write = csv.writer(csvfile)
			file_write.writerow(dist_list)
			csvfile.close()

		return 0

	except Exception as exception:
		print(exception)
		return -1

def read(name):
	try:
		f = open(name, 'r')
		csv_read = csv.reader(f, delimiter=';')
		matrix = []
		for line in csv_read:
			matrix.append(list(map(float, line[0].split(","))))
		return matrix
	except Exception as exception:
		print(exception)
		return -1

def read_sparse(name):
	try:
		f = open(name, 'r')
		csv_read = csv.reader(f, delimiter=';')
		row = []
		col = []
		data = []
		for line in csv_read:
			line_split = line[0].split(",")
			if len(line_split) == 3:
				row.append(int(line_split[0]))
				col.append(int(line_split[1]))
				data.append(float(line_split[2]))
			else:
				rows = int(line_split[0])
				cols = int(line_split[1])
		row = np.array(row)
		col = np.array(col)
		data = np.array(data)
		sparse_matrix = csr_matrix((data, (row, col)), shape=(rows, cols))
		return sparse_matrix

	except Exception as exception:
		print(exception)
		return -1

def parse_tree(tree_str):
	'''
	(Tint,lint,nodes_in_order) = parse_tree(tree_str)
	This function will parse a newick tree string and return the dictionary of ancestors Tint.
	Tint indexes the nodes by integers, Tint[i] = j means j is the ancestor of i.
	lint is a dictionary returning branch lengths: lint[(i,j)] = w(i,j) the weight of the edge connecting i and j.
	nodes_in_order is a list of the nodes in the input tree_str such that Tint[i]=j means nodes_in_order[j] is an ancestor
	of nodes_in_order[i]. Nodes are labeled from the leaves up.
	'''
	dtree = dendropy.Tree.get(data=tree_str, schema="newick", suppress_internal_node_taxa=False,store_tree_weights=True)
	#Name all the internal nodes
	nodes = dtree.nodes()
	i=0
	for node in nodes:
		if node.taxon == None:
			node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
			i = i+1
	full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]  # i in path from root to j only if i>j
	full_nodes_in_order.reverse()
	nodes_in_order = [item.taxon.label for item in full_nodes_in_order]  # i in path from root to j only if i>j
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = full_nodes_in_order[i]
		parent = node.parent_node
		if parent != None:
			Tint[i] = nodes_to_index[parent.taxon.label]
			lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
	return (Tint,lint,nodes_in_order)

def parse_tree_file(tree_str_file, suppress_internal_node_taxa=True, suppress_leaf_node_taxa=False):
	'''
	(Tint,lint,nodes_in_order) = parse_tree(tree_str_file)
	This function will parse a newick tree file (in the file given by tree_str_file) and return the dictionary of ancestors Tint.
	Tint indexes the nodes by integers, Tint[i] = j means j is the ancestor of i.
	lint is a dictionary returning branch lengths: lint[i,j] = w(i,j) the weight of the edge connecting i and j.
	nodes_in_order is a list of the nodes in the input tree_str such that T[i]=j means nodes_in_order[j] is an ancestor
	of nodes_in_order[i]. Nodes are labeled from the leaves up.
	'''
	dtree = dendropy.Tree.get(path=tree_str_file, schema="newick",
							suppress_internal_node_taxa=suppress_internal_node_taxa,
							store_tree_weights=True,
							suppress_leaf_node_taxa = suppress_leaf_node_taxa)
	#Name all the internal nodes
	nodes = dtree.nodes()
	i=0
	for node in nodes:
		if node.taxon == None:
			node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
			i = i+1
	full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]  # i in path from root to j only if i>j
	full_nodes_in_order.reverse()
	nodes_in_order = [item.taxon.label for item in full_nodes_in_order]  # i in path from root to j only if i>j
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = full_nodes_in_order[i]
		parent = node.parent_node
		if parent != None:
			Tint[i] = nodes_to_index[parent.taxon.label]
			if isinstance(node.edge.length, float):
				lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
			else:
				lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = 0.0
	return (Tint,lint,nodes_in_order)

def create_env(sample_file):
	'''
	:param sample_file: a file containding ids and samples
	:return: an env_dict in the form of { id: {sample:count} }
	'''
	env_dict = dict()
	with open(sample_file) as fp:
		line = fp.readline()
		while line:
			list = line.split()
			key = list.pop(0) #get key
			env_dict[key] = dict()
			for str in list:
				sample = str.split('_')[0]
				if sample in env_dict[key]:
					env_dict[key][sample] += 1
				else:
					env_dict[key][sample] = 1
			line = fp.readline()
	fp.close()
	return env_dict

def parse_envs(envs, nodes_in_order):
	'''
	(envs_prob_dict, samples) = parse_envs(envs, nodes_in_order)
	This function takes an environment envs and the list of nodes nodes_in_order and will return a dictionary envs_prob_dict
	with keys given by samples. envs_prob_dict[samples[i]] is a probability vector on the basis nodes_in_order denoting for sample i.
	'''
	nodes_in_order_dict = dict(zip(nodes_in_order,range(len(nodes_in_order))))
	for node in envs.keys():
		if node not in nodes_in_order_dict:
			print("Warning: environments contain taxa " + node + " not present in given taxonomic tree. Ignoring")
	envs_prob_dict = dict()
	for i in range(len(nodes_in_order)):
		node = nodes_in_order[i]
		if node in envs:
			samples = envs[node].keys()
			for sample in samples:
				if sample not in envs_prob_dict:
					envs_prob_dict[sample] = np.zeros(len(nodes_in_order))
					envs_prob_dict[sample][i] = envs[node][sample]
				else:
					envs_prob_dict[sample][i] = envs[node][sample]
	#Normalize samples
	samples = envs_prob_dict.keys()
	for sample in samples:
		if envs_prob_dict[sample].sum() == 0:
			warnings.warn("Warning: the sample %s has non-zero counts, do not use for Unifrac calculations" % sample)
		envs_prob_dict[sample] = envs_prob_dict[sample]/envs_prob_dict[sample].sum()
	return (envs_prob_dict, samples)