import sys, csv, itertools, os.path
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
from os import path
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics import rand_score
from sklearn_extra.cluster import KMedoids
import BiomWrapper as BW
import CSVWrapper as CSV
import MetadataWrapper as meta

if __name__ == "__main__":

	args = sys.argv
	if len(args) > 4:
		raise Exception("Invalid number of parameters.")
	elif len(args) == 4:
		L1_file = args[1]
		L2_file = args[2]
		num_clusters = int(args[3])
		if not path.exists(L1_file) or not path.exists(L2_file):
			raise Exception("Error: Invalid file path(s).")
	else:
		print("Error encountered in input arguments... Defaulting to L1-UniFrac-Out.csv and L2-UniFrac-Out.csv with 5 clusters... \n")
		L1_file = 'L1-UniFrac-Out.csv'
		L2_file = 'L2-UniFrac-Out.csv'
		num_clusters = 5
		if not path.exists(L1_file) or not path.exists(L2_file):
			raise Exception("Error: Missing default CSV file(s).")

	L1_distance_matrix = CSV.read(L1_file)

	f = open(L2_file, 'r')
	read = csv.reader(f, delimiter=';')
	L2_distance_matrix = []
	for i in read:
		L2_distance_matrix.append(list(map(float, i[0].split(","))))

	L1AgglomerativeCluster = AgglomerativeClustering(n_clusters=num_clusters, affinity='precomputed', linkage='complete').fit_predict(L1_distance_matrix)  
	L1KMedoidsCluster = KMedoids(n_clusters=num_clusters, metric='precomputed', method='pam', init='heuristic').fit_predict(L1_distance_matrix)

	L2AgglomerativeCluster = AgglomerativeClustering(n_clusters=num_clusters, affinity='precomputed', linkage='complete').fit_predict(L2_distance_matrix)  
	L2KMedoidsCluster = KMedoids(n_clusters=num_clusters, metric='precomputed', method='pam', init='heuristic').fit_predict(L2_distance_matrix)

	PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')
	metadata = meta.extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt')
	region_names = []
	for i in range(len(PCoA_Samples)):
		if metadata[PCoA_Samples[i]]['body_site'] not in region_names:
			region_names.append(metadata[PCoA_Samples[i]]['body_site'])
		PCoA_Samples[i] = region_names.index(metadata[PCoA_Samples[i]]['body_site'])

	print("Printing results for L1-UniFrac:")
	print("Metric\t\t\t\t\t\t\tAgglomerativeClustering\t\tKMedoids")
	RI1 = rand_score(PCoA_Samples, L1AgglomerativeCluster)
	RI2 = rand_score(PCoA_Samples, L1KMedoidsCluster)
	print(f"Rand Index Score:               {RI1}\t\t\t{RI2}")
	ARI1 = adjusted_rand_score(PCoA_Samples, L1AgglomerativeCluster)
	ARI2 = adjusted_rand_score(PCoA_Samples, L1KMedoidsCluster)
	print(f"Adjusted Rand Index Score:      {ARI1}\t\t\t{ARI2}")
	NMI1 = normalized_mutual_info_score(PCoA_Samples, L1AgglomerativeCluster)
	NMI2 = normalized_mutual_info_score(PCoA_Samples, L1KMedoidsCluster)
	print(f"Normalized Mutual Index Score:  {NMI1}\t\t\t{NMI2}")
	AMI1 = adjusted_mutual_info_score(PCoA_Samples, L1AgglomerativeCluster)
	AMI2 = adjusted_mutual_info_score(PCoA_Samples, L1KMedoidsCluster)
	print(f"Adjusted Mutual Info Score:     {AMI1}\t\t\t{AMI2}")
	FM1 = fowlkes_mallows_score(PCoA_Samples, L1AgglomerativeCluster)
	FM2 = fowlkes_mallows_score(PCoA_Samples, L1KMedoidsCluster)
	print(f"Fowlkes Mallows Score:          {FM1}\t\t\t{FM2}")

	print("\nPrinting results for L2-UniFrac:")
	print("Metric\t\t\t\t\t\t\tAgglomerativeClustering\t\tKMedoids")
	RI1 = rand_score(PCoA_Samples, L2AgglomerativeCluster)
	RI2 = rand_score(PCoA_Samples, L2KMedoidsCluster)
	print(f"Rand Index Score:               {RI1}\t\t\t{RI2}")
	ARI1 = adjusted_rand_score(PCoA_Samples, L2AgglomerativeCluster)
	ARI2 = adjusted_rand_score(PCoA_Samples, L2KMedoidsCluster)
	print(f"Adjusted Rand Index Score:      {ARI1}\t\t\t{ARI2}")
	NMI1 = normalized_mutual_info_score(PCoA_Samples, L2AgglomerativeCluster)
	NMI2 = normalized_mutual_info_score(PCoA_Samples, L2KMedoidsCluster)
	print(f"Normalized Mutual Index Score:  {NMI1}\t\t\t{NMI2}")
	AMI1 = adjusted_mutual_info_score(PCoA_Samples, L2AgglomerativeCluster)
	AMI2 = adjusted_mutual_info_score(PCoA_Samples, L2KMedoidsCluster)
	print(f"Adjusted Mutual Info Score:     {AMI1}\t\t\t{AMI2}")
	FM1 = fowlkes_mallows_score(PCoA_Samples, L2AgglomerativeCluster)
	FM2 = fowlkes_mallows_score(PCoA_Samples, L2KMedoidsCluster)
	print(f"Fowlkes Mallows Score:          {FM1}\t\t\t{FM2}")