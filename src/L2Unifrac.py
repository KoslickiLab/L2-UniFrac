import numpy as np
import dendropy
import sys
import warnings

epsilon = sys.float_info.epsilon

def parse_tree(tree_str):
	'''
	(Tint,lint,nodes_in_order) = parse_tree(tree_str)
	This function will parse a newick tree string and return the dictionary of ancestors Tint.
	Tint indexes the nodes by integers, Tint[i] = j means j is the ancestor of i.
	lint is a dictionary returning branch lengths: lint[(i,j)] = w(i,j) the weight of the edge connecting i and j.
	nodes_in_order is a list of the nodes in the input tree_str such that Tint[i]=j means nodes_in_order[j] is an ancestor
	of nodes_in_order[i]. Nodes are labeled from the leaves up.
	'''
	dtree = dendropy.Tree.get(data=tree_str, schema="newick", suppress_internal_node_taxa=False,store_tree_weights=True)
	#Name all the internal nodes
	nodes = dtree.nodes()
	i=0
	for node in nodes:
		if node.taxon == None:
			node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
			i = i+1
	full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]  # i in path from root to j only if i>j
	full_nodes_in_order.reverse()
	nodes_in_order = [item.taxon.label for item in full_nodes_in_order]  # i in path from root to j only if i>j
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = full_nodes_in_order[i]
		parent = node.parent_node
		if parent != None:
			Tint[i] = nodes_to_index[parent.taxon.label]
			lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
	return (Tint,lint,nodes_in_order)

def parse_tree_file(tree_str_file, suppress_internal_node_taxa=True, suppress_leaf_node_taxa=False):
	'''
	(Tint,lint,nodes_in_order) = parse_tree(tree_str_file)
	This function will parse a newick tree file (in the file given by tree_str_file) and return the dictionary of ancestors Tint.
	Tint indexes the nodes by integers, Tint[i] = j means j is the ancestor of i.
	lint is a dictionary returning branch lengths: lint[i,j] = w(i,j) the weight of the edge connecting i and j.
	nodes_in_order is a list of the nodes in the input tree_str such that T[i]=j means nodes_in_order[j] is an ancestor
	of nodes_in_order[i]. Nodes are labeled from the leaves up.
	'''
	dtree = dendropy.Tree.get(path=tree_str_file, schema="newick",
							suppress_internal_node_taxa=suppress_internal_node_taxa,
							store_tree_weights=True,
							suppress_leaf_node_taxa = suppress_leaf_node_taxa)
	#Name all the internal nodes
	nodes = dtree.nodes()
	i=0
	for node in nodes:
		if node.taxon == None:
			node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
			i = i+1
	full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]  # i in path from root to j only if i>j
	full_nodes_in_order.reverse()
	nodes_in_order = [item.taxon.label for item in full_nodes_in_order]  # i in path from root to j only if i>j
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = full_nodes_in_order[i]
		parent = node.parent_node
		if parent != None:
			Tint[i] = nodes_to_index[parent.taxon.label]
			if isinstance(node.edge.length, float):
				lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
			else:
				lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = 0.0
	return (Tint,lint,nodes_in_order)

def L2Unifrac_weighted_plain(ancestors, edge_lengths, nodes_in_order, P, Q):
	'''
	Z = EMDUnifrac_weighted(ancestors, edge_lengths, nodes_in_order, P, Q)
	This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
	and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
	Returns the weighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
	F[(i,j)] == num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node
	nodes_in_order[i] to the node nodes_in_order[j].
	'''
	num_nodes = len(nodes_in_order)
	Z = 0
	eps = 1e-8
	partial_sums = P - Q # Vector of partial sums obtained by computing the difference between probabilities of two samples. 
	#print(f"partial_sums: {partial_sums}")
	#print(f"ancestors: {ancestors}")
	total_mass = 1 

	for i in range(num_nodes - 1):
		print(f"iteration {i}'s Z: {Z}")
		val = partial_sums[i]
		print(f"partial_sums in loop: {partial_sums}")
		print(f"val: {val}")
		if abs(val) > eps:
			partial_sums[ancestors[i]] += val
			Z += edge_lengths[i, ancestors[i]]*(val**2)
	print(f"final Z: {Z}")
	Z = np.sqrt(Z)
	return Z

def push_up(P, Tint, lint, nodes_in_order):
	#print(P, Tint, lint, nodes_in_order)
	P_pushed = P + 0  # don't want to stomp on P
	for i in range(len(nodes_in_order) - 1):
		#print(i, Tint[i])
		if lint[i, Tint[i]] == 0:
			lint[i, Tint[i]] = epsilon
		#print(P_pushed)
		P_pushed[Tint[i]] += P_pushed[i]  # push mass up
		P_pushed[i] *= np.sqrt(lint[i, Tint[i]])
		#P_pushed[i] *= P_pushed[i]  # multiply mass at this node by edge length above it (TODO: P_push squared)
	return P_pushed

# 1) Push up and inverse push up on some vector and see if I get that vector back again.
# 2) Error most likely in either pushup or inverse pushup. Error is most likely in push up since test does not use inverse.
#    We could be potentially stomping on P_pushed as it propagates up.
#     a) Draw small tree, explicitly write mass and, by hand write what the pushup should be according to page 92/93 of thesis.
#	  b) Try permutations of push up to run test and see if it works.

def inverse_push_up(P, Tint, lint, nodes_in_order):
	P_pushed = np.zeros(P.shape)  # don't want to stomp on P
	for i in range(len(nodes_in_order) - 1):
		if lint[i, Tint[i]] == 0:
			edge_length = epsilon
		else:
			edge_length = lint[i, Tint[i]]
		p_val = P[i]
		P_pushed[i] += 1/np.sqrt(edge_length) * p_val  # re-adjust edge lengths
		if P_pushed[i] < epsilon:
			P_pushed[i] = 0
		P_pushed[Tint[i]] -= 1/np.sqrt(edge_length) * p_val  # propagate mass upward, via subtraction, only using immediate descendants
	root = len(nodes_in_order) - 1
	P_pushed[root] += P[root]
	return P_pushed

def mean_of_vectors(L):
	'''
	:param L: a list of vectors
	:return: a vector with each entry i being the mean of vectors of L at position i
	'''
	return np.mean(L, axis=0)

# COMPARE PUSHUP AND INVERSE PUSHUP WITH TESTS TO SEE THAT THEY CHANGE OR NOT

def create_env(sample_file):
	'''
	:param sample_file: a file containding ids and samples
	:return: an env_dict in the form of { id: {sample:count} }
	'''
	env_dict = dict()
	with open(sample_file) as fp:
		line = fp.readline()
		while line:
			list = line.split()
			key = list.pop(0) #get key
			env_dict[key] = dict()
			for str in list:
				sample = str.split('_')[0]
				if sample in env_dict[key]:
					env_dict[key][sample] += 1
				else:
					env_dict[key][sample] = 1
			line = fp.readline()
	fp.close()
	return env_dict

def parse_envs(envs, nodes_in_order):
	'''
	(envs_prob_dict, samples) = parse_envs(envs, nodes_in_order)
	This function takes an environment envs and the list of nodes nodes_in_order and will return a dictionary envs_prob_dict
	with keys given by samples. envs_prob_dict[samples[i]] is a probability vector on the basis nodes_in_order denoting for sample i.
	'''
	nodes_in_order_dict = dict(zip(nodes_in_order,range(len(nodes_in_order))))
	for node in envs.keys():
		if node not in nodes_in_order_dict:
			print("Warning: environments contain taxa " + node + " not present in given taxonomic tree. Ignoring")
	envs_prob_dict = dict()
	for i in range(len(nodes_in_order)):
		node = nodes_in_order[i]
		if node in envs:
			samples = envs[node].keys()
			for sample in samples:
				if sample not in envs_prob_dict:
					envs_prob_dict[sample] = np.zeros(len(nodes_in_order))
					envs_prob_dict[sample][i] = envs[node][sample]
				else:
					envs_prob_dict[sample][i] = envs[node][sample]
	#Normalize samples
	samples = envs_prob_dict.keys()
	for sample in samples:
		if envs_prob_dict[sample].sum() == 0:
			warnings.warn("Warning: the sample %s has non-zero counts, do not use for Unifrac calculations" % sample)
		envs_prob_dict[sample] = envs_prob_dict[sample]/envs_prob_dict[sample].sum()
	return (envs_prob_dict, samples)


def test_push_up():
	P1 = np.array([0.1, 0.2, 0,  0.3, 0, 0.3, 0.1])
	T1 = {0: 4, 1: 4, 2: 5, 3: 5, 4: 6, 5: 6}
	l1 = {(0, 4): 0.1, (1, 4): 0.1, (2, 5): 0.2, (3, 5): 0, (4, 6): 0.2, (5, 6): 0.2} # 0 edge_length not involving the root
	nodes1 = ['A', 'B', 'C', 'D', 'temp0', 'temp1', 'temp2']
	P_pushed1 = push_up(P1, T1, l1, nodes1)
	P_inversed1 = inverse_push_up(P_pushed1, T1, l1, nodes1)
	assert np.sum(abs(P1 - P_inversed1)) < 10**-10 #test inverse_push_up

#print(np.sum(abs(P1 - P_inversed1)))

#weight = L2Unifrac_weighted_plain(T1, l1, nodes1, P1, P1)
#print(weight)


def test_weighted_plain():
	tree_str = '((B:0.1,C:0.2)A:0.3);'  # there is an internal node (temp0) here.
	(T1, l1, nodes1) = parse_tree(tree_str)
	nodes_samples = {
	    'C': {'sample1': 1, 'sample2': 0},
	    'B': {'sample1': 1, 'sample2': 1},
	    'A': {'sample1': 0, 'sample2': 0},
	    'temp0': {'sample1': 0, 'sample2': 1}}  # temp0 is the root node
	(nodes_weighted, samples_temp) = parse_envs(nodes_samples, nodes1)
	#print(f"nodes weighted: {nodes_weighted}")
	#print(f"samples_temp: {samples_temp}")
	push_unifrac = np.linalg.norm(push_up(nodes_weighted['sample1'], T1, l1, nodes1) -
	              push_up(nodes_weighted['sample2'], T1, l1, nodes1))
	plain_unifrac = L2Unifrac_weighted_plain(T1, l1, nodes_samples, nodes_weighted['sample1'], nodes_weighted['sample2']) #calculated using L2Unifrac
	print(f"push_unifrac: {push_unifrac}")
	print(f"plain_unifrac: {plain_unifrac}")
	#print(f"plain_unifrac^2: {plain_unifrac**2}")
	#print(f"plain_unifrac sqrt: {np.sqrt(plain_unifrac)}")
	assert np.sum(np.abs(push_unifrac-plain_unifrac)) < 10**-10

def run_tests():
	test_push_up()
	test_weighted_plain()


if __name__ == '__main__':
	run_tests()