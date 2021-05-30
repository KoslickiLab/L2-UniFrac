import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L2Unifrac as L2U 
import BiomWrapper as BW
import dendropy
import numpy as np
import write_to_csv as CSV

if __name__ == "__main__":
	nodes_samples = BW.extract_biom('../data/47422_otu_table.biom')
	T1, l1, nodes_in_order = L2U.parse_tree_file('../data/trees/gg_13_5_otus_99_annotated.tree')
	(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)
	#print(samples_temp)
	L2UniFrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_in_order, nodes_weighted['1928.SRS058420.SRX020544.SRR045915'], nodes_weighted['1928.SRS011433.SRX020669.SRR045315'])
	print(L2UniFrac)

	PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')
	num_samples = len(PCoA_Samples)
	#print(PCoA_Samples)
	Distance_Matrix = []

	for i in range(num_samples):
		print(f"Iteration row: {i}")
		tmp_row = []
		for j in range(num_samples):
			if i < j:
				L2UniFrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_in_order, nodes_weighted[PCoA_Samples[i]], nodes_weighted[PCoA_Samples[j]])
				print(f"    Inner loop: {j} | L2-UniFrac: {L2UniFrac} | Sample 1: {PCoA_Samples[i]} | Sample 2: {PCoA_Samples[j]}")
			elif i > j:
				L2UniFrac = Distance_Matrix[j][i]
			else:
				L2UniFrac = 0.0
			tmp_row.append(L2UniFrac)
		
		Distance_Matrix.append(tmp_row)
		CSV.write('L2-UniFrac-Out.csv', tmp_row)

	#print(Distance_Matrix)