import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import CSVWrapper as CSV
import MetadataWrapper as meta
import numpy as np

PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')
metadata = meta.extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt')

region_names = []
region_map = {}
for i in range(len(PCoA_Samples)):
	if metadata[PCoA_Samples[i]]['body_site'] not in region_names:
		region_map[metadata[PCoA_Samples[i]]['body_site']] = []
		region_names.append(metadata[PCoA_Samples[i]]['body_site'])
	region_map[metadata[PCoA_Samples[i]]['body_site']].append(i)
	PCoA_Samples[i] = region_names.index(metadata[PCoA_Samples[i]]['body_site'])

sparse_matrix_L1 = CSV.read_sparse('L1-Push-Out.csv')

print(PCoA_Samples)
print(region_map)
for i in range(len(region_names)):
	group_arr = []
	for j in range(len(region_map[region_names[i]])):
		group_arr.append(np.array(sparse_matrix_L1[region_map[region_names[i]][j]].todense())[0])
	print(group_arr)
	print(len(group_arr))
	#average1 = L1U.median_of_vectors(sparse_matrix_L1.toarray())
	#print(average1)
	break

#sparse_matrix_L1 = CSV.read_sparse('L1-Push-Out.csv')
#sparse_matrix_L2 = CSV.read_sparse('L2-Push-Out.csv')

#average1 = L1U.median_of_vectors(sparse_matrix_L1.toarray())
#print(average1)