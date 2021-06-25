import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import csv, itertools
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score

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

	cluster = AgglomerativeClustering(n_clusters=num_clusters, affinity='precomputed', linkage='complete').fit_predict(distance_matrix)  

	#cluster = KMedoids(n_clusters=num_clusters, metric='precomputed', method='pam', init='heuristic').fit_predict(distance_matrix)

	print(cluster)

	cluster_arr = []
	for i in range(num_clusters):
		cluster_arr.append(np.where(cluster == i)[0])
	print(cluster_arr)

	PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')
	metadata = meta.extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt')
	region_names = []
	for i in range(len(PCoA_Samples)):
		if metadata[PCoA_Samples[i]]['body_site'] not in region_names:
			region_names.append(metadata[PCoA_Samples[i]]['body_site'])
		PCoA_Samples[i] = region_names.index(metadata[PCoA_Samples[i]]['body_site'])

	np.set_printoptions(threshold=sys.maxsize)
	print(cluster)
	print(PCoA_Samples)

	permutations = list(itertools.permutations(range(num_clusters), num_clusters))
	correct_max = 0
	correct_max_i = 0
	for i in range(len(permutations)):
		correct = 0
		set_map = {j:permutations[i][j] for j in range(num_clusters)}
		for j in range(len(cluster)):
			if set_map[cluster[j]] == PCoA_Samples[j]:
				correct += 1
		if correct > correct_max:
			correct_max = correct
			correct_max_i = i
	print(correct_max, correct_max/len(cluster), permutations[correct_max_i])
	# L1: 0.49217075984836
	# L2: 0.47651227954507996

	score = adjusted_rand_score(PCoA_Samples, cluster)
	print(score)
	score2 = normalized_mutual_info_score(PCoA_Samples, cluster)
	print(score2)
	score3 = adjusted_mutual_info_score(PCoA_Samples, cluster)
	print(score3)
	score4 = fowlkes_mallows_score(PCoA_Samples, cluster)
	print(score4)

	# L2:
	# 0.12190268522506367
	# 0.22313556847817295
	# 0.2222487106613287
	# 0.45534086122442907

	#count_dict = [[0 for j in range(num_clusters)] for i in range(num_clusters)]
	#for i in range(len(cluster)):
	#	count_dict[cluster[i]][PCoA_Samples[i]] += 1

	#disallowed_set = []
	#disallowed_match = []
	#translation_dict = {}
	#for iteration in range(num_clusters):
	#	print(count_dict)
	#	most_probable_set = None
	#	most_probable_match = None
	#	for i in range(num_clusters):
	#		for j in range(num_clusters):
	#			if most_probable_set == None and most_probable_match == None and i not in disallowed_set and j not in disallowed_match:
	#				most_probable_set = i
	#				most_probable_match = j
	#			elif i not in disallowed_set and j not in disallowed_match:
	#				if count_dict[i][j]/sum(count_dict[i]) > count_dict[most_probable_set][most_probable_match]/sum(count_dict[most_probable_set]) and i not in disallowed_set and j not in disallowed_match:
	#					most_probable_set = i
	#					most_probable_match = j
	#	print(count_dict[most_probable_set][most_probable_match]/sum(count_dict[most_probable_set]))
	#	disallowed_set.append(most_probable_set)
	#	disallowed_match.append(most_probable_match)
	#	translation_dict[most_probable_set] = most_probable_match
	#	print(most_probable_set, most_probable_match, disallowed_match) # 46.3%, 84.9%, 41.4%, 76.7%, 81.5%

	#num_correct = 0
	#num_incorrect = 0
	#for i in range(len(cluster)):
	#	cluster[i] = translation_dict[cluster[i]]
	#	if cluster[i] == PCoA_Samples[i]:
	#		num_correct += 1
	#	else:
	#		num_incorrect += 1

	#print(cluster)
	#print(num_correct, num_incorrect)
	#print(num_correct/len(cluster)) 

	# AgglomerativeClustering
	# 0.47008406131531233 - L1, 0.4720619746167793 - L2

	# KMedoids
	# 0.6245261249381903 - L1, 0.575243118509972 - L2

	#kmedoids = KMedoids(n_clusters=5, metric='precomputed', method='pam', init='heuristic').fit_predict(distance_matrix)
	#labels = kmedoids.labels_

	#print(kmedoids)