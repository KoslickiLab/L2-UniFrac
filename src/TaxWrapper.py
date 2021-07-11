def extract_tax(extension):
	try:
		with open(str(extension)) as f:
			taxonomies = {}
			tax_list = f.read().splitlines() 
			for tax in range(len(tax_list)):
				tax_list[tax] = tax_list[tax].split('\t')
				taxonomies[int(tax_list[tax][0])] = tax_list[tax][1]
			return taxonomies
	except:
		raise Exception("Unknown exception occurred. Try again.")

if __name__ == '__main__':
	print(extract_tax('../data/taxonomies/gg_13_8_99.gg.tax')[4545])\
	