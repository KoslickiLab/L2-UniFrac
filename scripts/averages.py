import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L1Unifrac as L1U
import L2Unifrac as L2U
import BiomWrapper as BW
import CSVWrapper as CSV
import MetadataWrapper as meta

sparse_matrix_L1 = CSV.read_sparse('L1-Push-Out.csv')
#sparse_matrix_L2 = CSV.read_sparse('L2-Push-Out.csv')

average1 = L1U.median_of_vectors(sparse_matrix_L1.toarray())
print(average1)