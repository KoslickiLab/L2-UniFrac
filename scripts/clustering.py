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
		file = 'L2-UniFrac-Out.csv'
	else:
		file = args[1]

	assert(isinstance(file, str))

	f = open(file, 'r')
	read = csv.reader(f, delimiter=';')
	distance_matrix = []
	for i in read:
		distance_matrix.append(list(map(float, i[0].split(","))))

	cluster = AgglomerativeClustering(n_clusters=5, affinity='precomputed', linkage = 'average').fit(distance_matrix)  

	print(cluster)

	print(np.where(cluster == 0)[0])
	print(np.where(cluster == 1)[0])
	print(np.where(cluster == 2)[0])
	print(np.where(cluster == 3)[0])
	print(np.where(cluster == 4)[0])

	print(np.where(cluster == 0)[0].size)
	print(np.where(cluster == 1)[0].size)
	print(np.where(cluster == 2)[0].size)
	print(np.where(cluster == 3)[0].size)
	print(np.where(cluster == 4)[0].size)

	#np.set_printoptions(threshold=sys.maxsize)

	#kmedoids = KMedoids(n_clusters=5, metric='precomputed', method='pam', init='heuristic').fit(distance_matrix)
	#labels = kmedoids.labels_