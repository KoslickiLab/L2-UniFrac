import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L2Unifrac as L2U 
import BiomWrapper as BW
import numpy as np

if __name__ == "__main__":
	nodes_samples = BW.extract_biom('../data/60982-reference-hit.biom')