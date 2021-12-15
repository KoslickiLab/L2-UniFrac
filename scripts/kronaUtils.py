import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import argparse, random, math, os, shutil, signal
from subprocess import STDOUT, TimeoutExpired, Popen, PIPE, run
from time import sleep

def generate_krona_visuals(region_names, tax_arr, inverse_pushed, output, intermediate_store=True):
	for name in region_names:
		tax_abundances = {}
		region_abundance_vector = inverse_pushed[name]
		for i in range(len(region_abundance_vector)):
			node_tax = tax_arr[i].split(';')
			new_node_tax = []
			for j in range(len(node_tax)):
				if node_tax[j] != 'root':
					new_node_tax.append(node_tax[j][3:])
				else:
					new_node_tax.append(node_tax[j])
			node_tax_str = '\t'.join(new_node_tax)
			if node_tax_str in tax_abundances:
				tax_abundances[node_tax_str] += region_abundance_vector[i]
			else:
				tax_abundances[node_tax_str] = region_abundance_vector[i]
		print(tax_abundances)
		for key in tax_abundances.keys():
			with open("krona/{0}_{1}_krona.txt".format(output, name), 'a') as file:
				file.write('\t'.join([str(tax_abundances[key]), key]))
		subpro = run('ktImportText krona/{0}_{1}_krona.txt -o krona/{0}_{1}.krona.html'.format(output, name), stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
		if not intermediate_store:
			os.remove('krona/{0}_{1}_krona.txt'.format(output, name))