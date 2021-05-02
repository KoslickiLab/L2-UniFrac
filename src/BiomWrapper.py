import biom
import numpy as np

def extract_biom(address):

	# Load table into the biom Table format.
	Biom = biom.load_table(address)

	# Grab the sample IDs (i.e. '1928.SRS045191.SRX021470.SRR052699')
	sample_ids = Biom.ids()

	print(sample_ids)

	# Grab the node IDs for the tree (i.e. '858026')
	phylogenetic_tree_nodes = Biom.ids(axis='observation')

	print(phylogenetic_tree_nodes)

	# Formulate the dictionary required for L2 Unifrac by associating each sample weight to its node ID.
	nodes_samples = {}
	for i in range(len(phylogenetic_tree_nodes)):
		row_vector = Biom._get_row(i).todense().tolist()[0] # Convert Table row to a standard python list.
		nodes_samples[phylogenetic_tree_nodes[i]] = {sample_ids[j]:row_vector[j] for j in range(len(sample_ids))}
	return nodes_samples

if __name__ == '__main__':
	nodes_sample = extract_biom('../data/60982-reference-hit.biom')
	#print(nodes_sample['4479984'])