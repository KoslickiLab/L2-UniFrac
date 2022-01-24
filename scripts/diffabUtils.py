import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import matplotlib.pyplot as plt
import L1Unifrac as L1U
import L2Unifrac as L2U

def generate_diffab(regions, region_averages, Tint, lint, nodes_in_order, taxonomy_in_order, output, thresh, maxDisp, includeTemp, L):
	if L == 1:
		for i in range(len(region_averages)):
			for j in range(len(region_averages)):
				if i < j:
					L1_UniFrac, DifferentialAbundance = L1U.EMDUnifrac_weighted(Tint, lint, nodes_in_order, region_averages[regions[i]], region_averages[regions[j]])
					fig = L2U.plot_diffab(nodes_in_order, taxonomy_in_order, DifferentialAbundance, regions[i], regions[j], plot_zeros=False, thresh=thresh, show=False, maxDisp=maxDisp, includeTemp=includeTemp)
					plt.savefig('images/{0}_diffab_{1}_{2}.png'.format(output, regions[i], regions[j]))
	else:
		for i in range(len(region_averages)):
			for j in range(len(region_averages)):
				if i < j:
					L2_UniFrac, DifferentialAbundance = L2U.L2Unifrac_weighted(Tint, lint, nodes_in_order, region_averages[regions[i]], region_averages[regions[j]])
					fig = L2U.plot_diffab(nodes_in_order, taxonomy_in_order, DifferentialAbundance, regions[i], regions[j], plot_zeros=False, thresh=thresh, show=False, maxDisp=maxDisp, includeTemp=includeTemp)
					plt.savefig('images/{0}_diffab_{1}_{2}.png'.format(output, regions[i], regions[j]))