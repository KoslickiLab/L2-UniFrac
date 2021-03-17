import numpy as np
import dendropy
import sys
import warnings

epsilon = sys.float_info.epsilon

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
	total_mass = 1 

	for i in range(num_nodes - 1):
		val = partial_sums[i]
		if abs(val) > eps:
			val2 = val**2
			partial_sums[ancestors[i]] += val
			Z += edge_lengths[i, ancestors[i]]*val2
	Z = np.sqrt(Z)
	return Z

def push_up(P, Tint, lint, nodes_in_order):
	P_pushed = P + 0  # don't want to stomp on P
	for i in range(len(nodes_in_order) - 1):
		if lint[i, Tint[i]] == 0:
			lint[i, Tint[i]] = epsilon
		P_pushed[Tint[i]] += P_pushed[i]  # push mass up
		P_pushed[i] *= lint[i, Tint[i]]  # multiply mass at this node by edge length above it (TODO: P_push squared)
	return P_pushed

def inverse_push_up(P, Tint, lint, nodes_in_order):
	P_pushed = np.zeros(P.shape)  # don't want to stomp on P
	for i in range(len(nodes_in_order) - 1):
		if lint[i, Tint[i]] == 0:
			edge_length = epsilon
		else:
			edge_length = lint[i, Tint[i]]
		p_val = P[i]
		P_pushed[i] += 1/edge_length * p_val  # re-adjust edge lengths (TODO: P_push squared/sqrt)
		P_pushed[Tint[i]] -= 1/edge_length * p_val  # propagate mass upward, via subtraction, only using immediate descendants
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