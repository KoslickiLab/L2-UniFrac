import numpy as np

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
