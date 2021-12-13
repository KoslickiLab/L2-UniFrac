import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import argparse, random, math, os, shutil, signal
from subprocess import STDOUT, TimeoutExpired, Popen, PIPE, run
from time import sleep

def generate_krona(region_names, tax_arr, inverse_pushed, output):
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
			with open("tmp_krona.txt", 'a') as file:
    			file.write('\t'.join([str(tax_abundances[key]), key]))
		subpro = run('ktImportText tmp_krona.txt -o {0}_{1}.krona.html'.format(output, name), stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
		os.remove('tmp_krona.txt')