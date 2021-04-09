import biom

def extract_biom(address):
	Biom = biom.load_table(address)
	print(Biom.ids()) 
	print(Biom.ids(axis='observation')) 
	print(Biom.nnz)
	#print(Biom)

if __name__ == '__main__':
	extract_biom('../data/47422_otu_table.biom')