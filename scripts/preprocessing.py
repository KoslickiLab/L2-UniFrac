import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import write_to_csv as CSV
import MetadataWrapper as meta
import multiprocessing as mp
import numpy

cores = mp.cpu_count()

nodes_samples = BW.extract_biom('../data/47422_otu_table.biom')
T1, l1, nodes_in_order = L2U.parse_tree_file('../data/trees/gg_13_5_otus_99_annotated.tree')
(nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes_in_order)

PCoA_Samples = BW.extract_samples('../data/47422_otu_table.biom')

import time
numpy.set_printoptions(threshold=sys.maxsize)
start = time.time()
L2_Pushed = L2U.push_up(nodes_weighted[PCoA_Samples[0]], T1, l1, nodes_in_order)
print(L2_Pushed)
end = time.time()
print(end-start)