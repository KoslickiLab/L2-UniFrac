import biom
import numpy as np

def extract_biom(address):
	Biom = biom.load_table(address)
	print(Biom.ids()) 
	phylogenetic_tree_nodes = Biom.ids(axis='observation')
	print(phylogenetic_tree_nodes) # Nodes in phylogenetic tree
	print(Biom.nnz)
	print(Biom)
	#print(Biom)
	#transform_f = lambda v,i,m: np.where(v % 3 == 0, v, 0)
	#mult_of_three = tform = Biom.transform(transform_f, inplace=False)
	#print(mult_of_three) 
	#phylum_idx = 1
	#collapse_f = lambda id_, md: '; '.join(md['taxonomy'][:phylum_idx + 1])
	#collapsed = mult_of_three.collapse(collapse_f, axis='observation')
	#print(collapsed) 
	#pa = collapsed.pa()
	#print(pa)
	#print(Biom)

if __name__ == '__main__':
	extract_biom('../data/sampleBiom.biom')