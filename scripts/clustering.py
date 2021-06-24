import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import csv
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn_extra.cluster import KMedoids
import BiomWrapper as BW
import metadata_wrapper as meta
import matplotlib.pyplot as plt

if __name__ == "__main__":

	args = sys.argv
	if len(args) > 3:
		raise Exception("Invalid number of parameters.")
	elif len(args) == 1:
		file = 'L1-UniFrac-Out.csv'
		num_clusters = 5
	elif len(args) == 2:
		file = args[1]
		num_clusters = 5
	else:
		file = args[1]
		num_clusters = int(args[2])

	f = open(file, 'r')
	read = csv.reader(f, delimiter=';')
	distance_matrix = []
	for i in read:
		distance_matrix.append(list(map(float, i[0].split(","))))

	cluster = AgglomerativeClustering(n_clusters=num_clusters, affinity='precomputed', linkage = 'complete').fit_predict(distance_matrix)  

	print(cluster)

	cluster_arr = []
	for i in range(num_clusters):
		cluster_arr.append(np.where(cluster == i)[0])
	print(cluster_arr)

	PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')
	metadata = meta.extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt')
	for i in range(len(PCoA_Samples)):
		PCoA_Samples[i] = metadata[PCoA_Samples[i]]['body_site']

	print(PCoA_Samples)

	#for i in range(len(cluster_arr)):
		#for j in range(len(cluster_arr[i])):


	#np.set_printoptions(threshold=sys.maxsize)

	#kmedoids = KMedoids(n_clusters=5, metric='precomputed', method='pam', init='heuristic').fit(distance_matrix)
	#labels = kmedoids.labels_