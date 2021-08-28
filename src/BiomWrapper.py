import biom
import numpy as np

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

if __name__ == '__main__':
	
	# Test biom
	nodes_test = extract_biom('../data/47422_otu_table.biom')
	print(len(nodes_test.keys()))
	print(nodes_test['1050608'])
	print(sum(nodes_test['1050608'].values()))
	tot = 0
	Biom = biom.load_table('../data/47422_otu_table.biom')
	sample_ids = Biom.ids()
	phylogenetic_tree_nodes = Biom.ids(axis='observation')
	for i in range(len(phylogenetic_tree_nodes)):
		row = nodes_test[phylogenetic_tree_nodes[i]]
		for index, (key, value) in enumerate(row.items()):
			tot += nodes_test[phylogenetic_tree_nodes[i]][key]
			break
	print("Total: ", tot)
	#print(len(extract_samples('../data/47422_otu_table.biom')))

	# Real-world biom
	#nodes_sample = extract_biom('../data/47422_otu_table.biom')
	#print(nodes_sample[:8])
	#print(len(nodes_sample))