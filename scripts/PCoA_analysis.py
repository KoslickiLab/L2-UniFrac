import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
from skbio.stats.ordination import pcoa
import csv
import BiomWrapper as BW
import pandas as pd
from skbio import DistanceMatrix
import matplotlib.pyplot as plt

f = open('L2-UniFrac-Out.csv', 'r')
read = csv.reader(f, delimiter=';')
distance_matrix = []
for i in read:
	distance_matrix.append(list(map(float, i[0].split(","))))

sk_distance_matrix = DistanceMatrix(distance_matrix, BW.extract_samples('../data/47422_otu_table.biom')[:8])

metadata = {'1928.SRS015139.SRX020561.SRR043887': {'site': '1'},
			'1928.SRS015198.SRX020561.SRR043886': {'site': '2'},
			'1928.SRS015198.SRX020582.SRR049873': {'site': '3'},
			'1928.SRS019552.SRX020561.SRR043896': {'site': '4'},
			'1928.SRS019556.SRX020561.SRR043882': {'site': '5'},
			'1928.SRS019558.SRX020561.SRR043890': {'site': '6'},
			'1928.SRS019679.SRX020561.SRR043905': {'site': '7'},
			'1928.SRS046325.SRX022229.SRR058092': {'site': '8'},
			}

pd_metadata = pd.DataFrame.from_dict(metadata, orient='index')

result = pcoa(sk_distance_matrix)
print(result)

fig = result.plot(df=pd_metadata, column='site',
						title='Samples colored by body site',
						cmap='Set1', s=50)

plt.show()