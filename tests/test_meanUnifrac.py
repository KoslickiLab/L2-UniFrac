import sys
sys.path.append('../L2Unifrac')
sys.path.append('../L2Unifrac/src')
sys.path.append('../src')
import L2Unifrac as L2U 
import numpy as np

try:
    (Tint, lint, nodes_in_order) = L2U.parse_tree_file('../data/old_UniFrac/97_otus_unannotated.tree')
    env_dict = L2U.create_env('../data/old_UniFrac/289_seqs_otus.txt')
except FileNotFoundError:
    (Tint, lint, nodes_in_order) = L2U.parse_tree_file('../data/old_UniFrac/97_otus_unannotated.tree')
    env_dict = L2U.create_env('../data/old_UniFrac/289_seqs_otus.txt')
(env_prob_dict, samples) = L2U.parse_envs(env_dict, nodes_in_order)

#test parse_tree
def test_parse_tree():
    tree_str = '((B:0.1,C:0.2)A:0.3);'
    (Tint1, lint1, nodes_in_order1) = L2U.parse_tree(tree_str)
    assert Tint1 == {0: 2, 1: 2, 2: 3}
    assert lint1 == {(1, 2): 0.1, (2, 3): 0.3, (0, 2): 0.2}
    assert nodes_in_order1 == ['C', 'B', 'A', 'temp0']  # temp0 is the root node

#test push_up and inverse_push_up
def test_inverse():
    #simple tests
    P1 = np.array([0.1, 0.2, 0,  0.3, 0, 0.3, 0.1])
    T1 = {0: 4, 1: 4, 2: 5, 3: 5, 4: 6, 5: 6}
    l1 = {(0, 4): 0.1, (1, 4): 0.1, (2, 5): 0.2, (3, 5): 0, (4, 6): 0.2, (5, 6): 0.2} # 0 edge_length not involving the root
    nodes1 = ['A', 'B', 'C', 'D', 'temp0', 'temp1', 'temp2']
    P_pushed1 = L2U.push_up(P1, T1, l1, nodes1)
    x = np.sqrt(L2U.epsilon) * 0.3
    answer1 = np.array([0.0316227766, 0.0632455532, 0, x, 0.134164079, 0.268328157, 1])
    assert all(np.abs(P_pushed1 - answer1) < 0.00000001) #test push_up
    assert P_pushed1[3] > 10**-18 #P_pushed[3] (edge length 0) is non-zero
    P_inversed1 = L2U.inverse_push_up(P_pushed1, T1, l1, nodes1)
    assert np.sum(abs(P1 - P_inversed1)) < 10**-10 #test inverse_push_up
    l2 = {(0, 4): 0.1, (1, 4): 0.1, (2, 5): 0.2, (3, 5): 0, (4, 6): 0.2, (5, 6): 0} # more than one edge with 0 edge length, involving the root.
    y = np.sqrt(L2U.epsilon) * 0.6
    P_pushed2 = L2U.push_up(P1, T1, l2, nodes1)
    answer2 = np.array([0.0316227766, 0.0632455532, 0, x, 0.134164079, y, 1])
    assert all(np.abs(P_pushed2 - answer2) < 0.00000001)
    #test with real data
    Q = env_prob_dict['232.M2Lsft217']
    Q_pushed = L2U.push_up(Q, Tint, lint, nodes_in_order)
    Q_inversed = L2U.inverse_push_up(Q_pushed, Tint, lint, nodes_in_order)
    assert np.sum(abs(Q - Q_inversed)) < 10**-10

#test if push_up computes the correct unifrac value
def test_push_up():
    tree_str = '((B:0.1,C:0.2)A:0.3);'  # there is an internal node (temp0) here.
    (T1, l1, nodes1) = L2U.parse_tree(tree_str)
    nodes_samples = {
        'C': {'sample1': 1, 'sample2': 0},
        'B': {'sample1': 1, 'sample2': 1},
        'A': {'sample1': 0, 'sample2': 0},
        'temp0': {'sample1': 0, 'sample2': 1}}  # temp0 is the root node
    (nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes1)
    unifrac2 = np.linalg.norm(L2U.push_up(nodes_weighted['sample1'], T1, l1, nodes1) -
                  L2U.push_up(nodes_weighted['sample2'], T1, l1, nodes1))
    #EMDUnifrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_samples, nodes_weighted['sample1'], nodes_weighted['sample2']) #calculated using L2Unifrac
    #print(unifrac2, EMDUnifrac)
    #assert np.abs(unifrac2 - EMDUnifrac) < 10**-8

    tree_str = '((B:0.1,C:0.2)A:0.3);'  # there is an internal node (temp0) here.
    (T1, l1, nodes1) = L2U.parse_tree(tree_str)
    nodes_samples = {
        'C': {'sample1': 1, 'sample2': 0},
        'B': {'sample1': 1, 'sample2': 1},
        'A': {'sample1': 1, 'sample2': 0},
        'temp0': {'sample1': 0, 'sample2': 1}}  # temp0 is the root node
    (nodes_weighted, samples_temp) = L2U.parse_envs(nodes_samples, nodes1)
    print(nodes_weighted)
    unifrac2 = np.linalg.norm(L2U.push_up(nodes_weighted['sample1'], T1, l1, nodes1) -
                  L2U.push_up(nodes_weighted['sample2'], T1, l1, nodes1))
    EMDUnifrac = L2U.L2Unifrac_weighted_plain(T1, l1, nodes_samples, nodes_weighted['sample1'], nodes_weighted['sample2']) #calculated using L2Unifrac
    print(unifrac2, EMDUnifrac)
    assert np.abs(unifrac2 - EMDUnifrac) < 10**-8
    #assert np.sum(np.abs(unifrac2 - EMDUnifrac)) < 10**-10
    #assert unifrac1 == 0.25
    #test with real data
    P = env_prob_dict['232.M9Okey217']
    Q = env_prob_dict['232.M3Indl217']
    unifrac2 = np.linalg.norm(L2U.push_up(P, Tint, lint, nodes_in_order) -
                             L2U.push_up(Q, Tint, lint, nodes_in_order))
    EMDUnifrac2 = L2U.L2Unifrac_weighted_plain(Tint, lint, env_dict, P, Q) #calculated using L2Unifrac
    print(unifrac2, EMDUnifrac2)
    assert np.abs(unifrac2 - EMDUnifrac2) < 10**-8
    #assert np.sum(np.abs(unifrac2 - EMDUnifrac)) < 10**-10

def run_tests():
    test_parse_tree()
    test_inverse()
    test_push_up()

run_tests()