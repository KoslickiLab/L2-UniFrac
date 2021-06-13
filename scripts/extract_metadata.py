import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import BiomWrapper as BW
from itertools import islice

# Function to pull metadata from specific metadata file sourced at http://mse.ac.cn/index.php/mse
def extract_metadata(extension):
	try:
		f = open(str(extension), 'r')
		metadata = f.readlines()

		for i in range(len(metadata)):
			metadata[i] = metadata[i].split("\t")

		del[metadata[0]]

		meta_dict = {}
		for i in range(len(metadata)):
			meta_dict[metadata[i][0]] = metadata[i][3][7:]

		return meta_dict

	except FileNotFoundError:
		raise FileNotFoundError("Unknown file. Make sure you have the correct address and try again!")

	except:
		raise Exception("Unknown exception occurred. Try again.")

# Ensures that extracted metadata is suffient for the selected biom file
def test_metadata_completeness(metadata_path, biom_extension):
	meta_dict = extract_metadata(metadata_path)
	samples = BW.extract_samples(biom_extension)

	for i in range(len(samples)):
		assert samples[i] in meta_dict

if __name__ == "__main__":

	# Ensures that the metadata corresponds with the sample data selection
	test_metadata_completeness('../data/metadata/P_1928_65684500_raw_meta.txt', '../data/47422_otu_table.biom')

	# 4 items from the dictionary as an example
	print(list(islice(extract_metadata('../data/metadata/P_1928_65684500_raw_meta.txt').items(), 4)))