import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
from os import path
from skbio.stats.ordination import pcoa
import csv
import BiomWrapper as BW
import CSVWrapper as CSV
import MetadataWrapper as meta
import pandas as pd
from skbio import DistanceMatrix
import matplotlib.pyplot as plt

# File cheatsheet total (python PCoA_analysis.py 0 L2-UniFrac-Out.csv ../data/47422_otu_table.biom ../data/metadata/P_1928_65684500_raw_meta.txt):
# option:		 
# distance_file: 'L2-UniFrac-Out.csv'
# biom_file:     '../data/47422_otu_table.biom'
# metadata_file: '../data/metadata/P_1928_65684500_raw_meta.txt'

# File cheatsheet groups (python PCoA_analysis.py 1 L2-Group-UniFrac-Out.csv ../data/47422_otu_table.biom ../data/metadata/P_1928_65684500_raw_meta.txt skin,saliva,oral\ cavity,vagina,feces):
# option:		 
# distance_file: 'L2-UniFrac-Out.csv'
# biom_file:     '../data/47422_otu_table.biom'
# metadata_file: '../data/metadata/P_1928_65684500_raw_meta.txt'

def PCoA_total(distance_file, biom_file, metadata_file, plot=True):
	distance_matrix = CSV.read(distance_file)

	sk_distance_matrix = DistanceMatrix(distance_matrix, BW.extract_samples(biom_file))

	metadata = meta.extract_metadata(metadata_file)

	pd_metadata = pd.DataFrame.from_dict(metadata, orient='index')

	result = pcoa(sk_distance_matrix)
	print(result)

	fig = result.plot(df=pd_metadata, column='body_site',
							axis_labels=('PC 1', 'PC 2', 'PC 3'),
							title='Samples colored by body site',
							cmap='Set1', s=50)

	if plot:
		plt.show()
	else:
		return fig

def PCoA_group(distance_file, biom_file, groups, plot=True):
	distance_matrix = CSV.read(distance_file)

	sk_distance_matrix = DistanceMatrix(distance_matrix, [str(i) for i in range(len(groups))])

	metadata = {str(i): {'body_site': groups[i]} for i in range(len(groups))}

	pd_metadata = pd.DataFrame.from_dict(metadata, orient='index')

	result = pcoa(sk_distance_matrix)
	print(result)

	fig = result.plot(df=pd_metadata, column='body_site',
							axis_labels=('PC 1', 'PC 2', 'PC 3'),
							title='Samples colored by body site',
							cmap='Set1', s=50)

	if plot:
		plt.show()
	else:
		return fig

def PCoA_total_from_matrix(distance_matrix, biom_file, metadata_file, plot=True):
	sk_distance_matrix = DistanceMatrix(distance_matrix, BW.extract_samples(biom_file))

	metadata = meta.extract_metadata(metadata_file)

	pd_metadata = pd.DataFrame.from_dict(metadata, orient='index')

	result = pcoa(sk_distance_matrix)
	print(result)

	fig = result.plot(df=pd_metadata, column='body_site',
							axis_labels=('PC 1', 'PC 2', 'PC 3'),
							title='Samples colored by body site',
							cmap='Set1', s=50)

	if plot:
		plt.show()
	else:
		return fig

def PCoA_group_from_matrix(distance_matrix, biom_file, groups, plot=True):
	sk_distance_matrix = DistanceMatrix(distance_matrix, [str(i) for i in range(len(groups))])

	metadata = {str(i): {'body_site': groups[i]} for i in range(len(groups))}

	pd_metadata = pd.DataFrame.from_dict(metadata, orient='index')

	result = pcoa(sk_distance_matrix)
	print(result)

	fig = result.plot(df=pd_metadata, column='body_site',
							axis_labels=('PC 1', 'PC 2', 'PC 3'),
							title='Samples colored by body site',
							cmap='Set1', s=50)

	if plot:
		plt.show()
	else:
		return fig

if __name__ == "__main__":
	args = sys.argv
	if len(args) != 5 and len(args) != 6:
		raise Exception("Invalid number of parameters.")
	else:
		option = args[1]
		distance_file = args[2]
		biom_file = args[3]
		metadata_file = args[4]
		print(distance_file, biom_file, metadata_file)
		if not path.exists(distance_file) or not path.exists(biom_file) or not path.exists(metadata_file):
			raise Exception("Error: Invalid file path(s).")

		if int(option) == 0:
			PCoA_total(distance_file, biom_file, metadata_file)
		elif int(option) == 1:
			body_sites = args[5].split(",")
			PCoA_group(distance_file, biom_file, body_sites)