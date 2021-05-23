import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L2Unifrac as L2U 
import BiomWrapper as BW
import dendropy
import numpy as np

if __name__ == "__main__":
	nodes_samples = BW.extract_biom('../data/47422_otu_table.biom')
	T1, l1, nodes1 = L2U.parse_tree_file('../data/trees/gg_99_otus_4feb2011.tre')
	(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes1)
	EMDUnifrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_samples, nodes_weighted['1928.SRS058420.SRX020544.SRR045915'], nodes_weighted['1928.SRS011433.SRX020669.SRR045315'])
	print(EMDUnifrac)