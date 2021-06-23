import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import csv
import numpy as np
from sklearn.cluster import AgglomerativeClustering

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

	#sk_distance_matrix = DistanceMatrix(distance_matrix, BW.extract_samples('../data/47422_otu_table.biom'))

	agg = AgglomerativeClustering(n_clusters=5, affinity='precomputed', linkage = 'complete')
	agg2 = agg.fit_predict(distance_matrix)  # Returns class labels.

	print(agg)
	print(agg2)

	print(np.where(agg2 == 0))
	print(np.where(agg2 == 1))
	print(np.where(agg2 == 2))
	print(np.where(agg2 == 3))
	print(np.where(agg2 == 4))

	#metadata = meta.extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt')

	#pd_metadata = pd.DataFrame.from_dict(metadata, orient='index')

	#result = pcoa(sk_distance_matrix)
	#print(result)

	#fig = result.plot(df=pd_metadata, column='body_site',
	#						axis_labels=('PC 1', 'PC 2', 'PC 3'),
	#						title='Samples colored by body site',
	#						cmap='Set1', s=50)

	#plt.show()