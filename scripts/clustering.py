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

def report_clustering(distance_file, biom_file, metadata_file, verbose, L=2, output_file=None):
	if not isinstance(distance_file, list):
		distance_matrix = CSV.read(distance_file)
	else:
		distance_matrix = distance_file

	if output_file is not None:
		f = open(output_file, 'w')

	output_matrix = []

	AgglomerativeCluster = AgglomerativeClustering(n_clusters=num_clusters, affinity='precomputed', linkage='complete').fit_predict(distance_matrix)  
	KMedoidsCluster = KMedoids(n_clusters=num_clusters, metric='precomputed', method='pam', init='heuristic').fit_predict(distance_matrix)

	PCoA_Samples = BW.extract_samples(biom_file)
	metadata = meta.extract_metadata(metadata_file)
	region_names = []
	for i in range(len(PCoA_Samples)):
		if metadata[PCoA_Samples[i]]['body_site'] not in region_names:
			region_names.append(metadata[PCoA_Samples[i]]['body_site'])
		PCoA_Samples[i] = region_names.index(metadata[PCoA_Samples[i]]['body_site'])

	if verbose and L == 1:
		print('Printing results for L1-UniFrac:')
	elif verbose and L == 2:
		print('Printing results for L2-UniFrac:')
	if verbose:
		print('Metric\t\t\t\t\t\t\tAgglomerativeClustering\t\tKMedoids')

	if output_file is not None:
		if L == 1:
			f.write('Printing results for L1-UniFrac:\n')
		elif L == 2:
			f.write('Printing results for L2-UniFrac:\n')
		f.write('Metric\t\t\t\tAgglomerativeClustering\t\t\tKMedoids\n')

	if L == 1:
		output_matrix.append(['Printing results for L1-UniFrac:'])
	if L == 2:
		output_matrix.append(['Printing results for L2-UniFrac:'])
	output_matrix.append(['Metric', 'AgglomerativeClustering', 'KMedoids'])

	RI1 = rand_score(PCoA_Samples, AgglomerativeCluster)
	RI2 = rand_score(PCoA_Samples, KMedoidsCluster)
	if verbose:
		print(f'Rand Index Score:               {RI1}\t\t\t{RI2}')
	ARI1 = adjusted_rand_score(PCoA_Samples, AgglomerativeCluster)
	ARI2 = adjusted_rand_score(PCoA_Samples, KMedoidsCluster)
	if verbose:
		print(f'Adjusted Rand Index Score:      {ARI1}\t\t\t{ARI2}')
	NMI1 = normalized_mutual_info_score(PCoA_Samples, AgglomerativeCluster)
	NMI2 = normalized_mutual_info_score(PCoA_Samples, KMedoidsCluster)
	if verbose:
		print(f'Normalized Mutual Index Score:  {NMI1}\t\t\t{NMI2}')	
	AMI1 = adjusted_mutual_info_score(PCoA_Samples, AgglomerativeCluster)
	AMI2 = adjusted_mutual_info_score(PCoA_Samples, KMedoidsCluster)
	if verbose:
		print(f'Adjusted Mutual Info Score:     {AMI1}\t\t\t{AMI2}')
	FM1 = fowlkes_mallows_score(PCoA_Samples, AgglomerativeCluster)
	FM2 = fowlkes_mallows_score(PCoA_Samples, KMedoidsCluster)
	if verbose:
		print(f'Fowlkes Mallows Score:          {FM1}\t\t\t{FM2}')

	if output_file is not None:
		f.write(f'Rand Index Score:               {RI1}\t\t\t{RI2}\n')
		f.write(f'Adjusted Rand Index Score:      {ARI1}\t\t\t{ARI2}\n')
		f.write(f'Normalized Mutual Index Score:  {NMI1}\t\t\t{NMI2}\n')
		f.write(f'Adjusted Mutual Info Score:     {AMI1}\t\t\t{AMI2}\n')
		f.write(f'Fowlkes Mallows Score:          {FM1}\t\t\t{FM2}\n')

	output_matrix.append(['Rand Index Score:', RI1, RI2])
	output_matrix.append(['Adjusted Rand Index Score:', ARI1, ARI2])
	output_matrix.append(['Normalized Mutual Index Score:', NMI1, NMI2])
	output_matrix.append(['Adjusted Mutual Info Score:', AMI1, AMI2])
	output_matrix.append(['Fowlkes Mallows Score:', FM1, FM2])

	return output_matrix

if __name__ == '__main__':

	args = sys.argv
	if len(args) > 4:
		raise Exception('Invalid number of parameters.')
	elif len(args) == 4:
		L1_file = args[1]
		L2_file = args[2]
		num_clusters = int(args[3])
		if not path.exists(L1_file) or not path.exists(L2_file):
			raise Exception('Error: Invalid file path(s).')
	else:
		print('Error encountered in input arguments... Defaulting to L1-UniFrac-Out.csv and L2-UniFrac-Out.csv with 5 clusters... \n')
		L1_file = 'L1-UniFrac-Out.csv'
		L2_file = 'L2-UniFrac-Out.csv'
		num_clusters = 5
		if not path.exists(L1_file) or not path.exists(L2_file):
			raise Exception('Error: Missing default CSV file(s).')

	print(report_clustering(L1_file, '../data/47422_otu_table.biom', '../data/metadata/P_1928_65684500_raw_meta.txt', False, 1))
	print(report_clustering(L2_file, '../data/47422_otu_table.biom', '../data/metadata/P_1928_65684500_raw_meta.txt', False, 2))
