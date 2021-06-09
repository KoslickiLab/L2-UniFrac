import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L2Unifrac as L2U
import BiomWrapper as BW
import write_to_csv as CSV
import multiprocessing

cores = multiprocessing.cpu_count()

nodes_samples = BW.extract_biom('../data/47422_otu_table.biom')
T1, l1, nodes_in_order = L2U.parse_tree_file('../data/trees/gg_13_5_otus_99_annotated.tree')
(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)
#print(samples_temp)
L2UniFrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_in_order, nodes_weighted['1928.SRS058420.SRX020544.SRR045915'], nodes_weighted['1928.SRS011433.SRX020669.SRR045315'])
print(L2UniFrac)

PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')

Distance_Matrix = []

# Testing subset of samples...
PCoA_Samples = PCoA_Samples[:512]

def unifrac_work_wrapper(args):
	#print(args)
	return unifrac_worker(*args)

def unifrac_worker(samp1num, samp2num):
	#print(f"Samples: {samp1num}, {samp2num}")
	#print(f"PCoA_Samples: {PCoA_Samples}")
	L2UniFrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_in_order, nodes_weighted[PCoA_Samples[samp1num]], nodes_weighted[PCoA_Samples[samp2num]])
	formatted_L2 = "{:.16f}".format(L2UniFrac)
	return L2UniFrac, f"\tInner loop: {str(samp2num).zfill(4)} | L2-UniFrac: {formatted_L2} | Sample 1: {PCoA_Samples[samp1num]} | Sample 2: {PCoA_Samples[samp2num]}"

if __name__ == "__main__":

	args = sys.argv
	if len(args) > 1:
		if args[1] == "1":
			debug = 1
		else:
			debug = 0
	else:
		debug = 0

	if debug == 1:
		local_vars = list(locals().items())
		for var, obj in local_vars:
			print(f"{var.ljust(17)}: {sys.getsizeof(obj)}")

	# Multi Core Method
	#for i in range(len(PCoA_Samples)):

	#        print(f"Iteration row: {i}")

	row = [(i, j) for j in range(len(PCoA_Samples)) for i in range(len(PCoA_Samples))]

	with multiprocessing.Pool(processes=cores-1) as pool:
		result = pool.map(unifrac_work_wrapper, row)

	for i in range(len(PCoA_Samples)):
		dist_list = []
		for j in range(len(PCoA_Samples)):
			dist_list.append(result[i*len(PCoA_Samples)+j][0])
			if debug == 1:
				print(result[i*len(PCoA_Samples)+j][1])

		CSV.write('L2-UniFrac-Out.csv', dist_list)

	# Single Core Method
	#for i in range(len(PCoA_Samples)):
	#	print(f"Iteration row: {i}")
	#	tmp_row = []
	#	for j in range(len(PCoA_Samples)):
	#		#if i < j:
	#		L2UniFrac, out = unifrac_worker(i, j)
	#		#print(f"    Inner loop: {j} | L2-UniFrac: {L2UniFrac} | Sample 1: {PCoA_Samples[i]} | Sample 2: {PCoA_Samples[j]}")
	#		print(out)
	#		#elif i > j:
	#		#	L2UniFrac = Distance_Matrix[j][i]
	#		#else:
	#		#	L2UniFrac = 0.0
	#		tmp_row.append(L2UniFrac)
	#	
	#	Distance_Matrix.append(tmp_row)
	#	CSV.write('L2-UniFrac-Out.csv', tmp_row)

	#print(Distance_Matrix)