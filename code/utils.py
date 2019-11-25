#!/usr/bin/env python

# Utility functions

import logging
import sys
import networkx as nx
import numpy as np
from scipy import sparse
from scipy.linalg import eig
import matplotlib.pyplot as plt
from itertools import product


'''
Functions get_base_modularity_matrix, get_modularity, get_mod_matrix, largest_eig from https://github.com/zhiyzuo/python-modularity-maximization


License notice:

Copyright 2018 Zhiya Zuo

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

'''

def get_base_modularity_matrix(network):
    '''
    Obtain the modularity matrix for the whole network

    Parameters
    ----------
    network : nx.Graph or nx.DiGraph
        The network of interest

    Returns
    -------
    np.matrix
        The modularity matrix for `network`

    Raises
    ------
    TypeError
        When the input `network` does not fit either nx.Graph or nx.DiGraph
    '''

    if type(network) == nx.Graph:
        #return sparse.csc_matrix(nx.modularity_matrix(network))
        return nx.modularity_matrix(network)
    else:
        raise TypeError('Graph type not supported. Use nx.Graph')

def get_modularity(network, community_dict):
    '''
    Calculate the modularity. Edge weights are ignored.

    Undirected:
    .. math:: Q = \frac{1}{2m}\sum_{i,j} \(A_ij - \frac{k_i k_j}{2m}\) * \detal_(c_i, c_j)

    Directed:
    .. math:: Q = \frac{1}{m}\sum_{i,j} \(A_ij - \frac{k_i^{in} k_j^{out}}{m}\) * \detal_{c_i, c_j}

    Parameters
    ----------
    network : nx.Graph or nx.DiGraph
        The network of interest
    community_dict : dict
        A dictionary to store the membership of each node
        Key is node and value is community index

    Returns
    -------
    float
        The modularity of `network` given `community_dict`
    '''

    Q = 0
    G = network.copy()
    nx.set_edge_attributes(G, {e:1 for e in G.edges}, 'weight')
    A = nx.to_scipy_sparse_matrix(G).astype(float)

    if type(G) == nx.Graph:
        # for undirected graphs, in and out treated as the same thing
        out_degree = in_degree = dict(nx.degree(G))
        M = 2.*(G.number_of_edges())
        print("Calculating modularity for undirected graph")
    elif type(G) == nx.DiGraph:
        in_degree = dict(G.in_degree())
        out_degree = dict(G.out_degree())
        M = 1.*G.number_of_edges()
        print("Calculating modularity for directed graph")
    else:
        print('Invalid graph type')
        raise TypeError

    nodes = list(G)
    Q = np.sum([A[i,j] - in_degree[nodes[i]]*\
                         out_degree[nodes[j]]/M\
                 for i, j in product(range(len(nodes)),\
                                     range(len(nodes))) \
                if community_dict[nodes[i]] == community_dict[nodes[j]]])
    return Q / M


def get_mod_matrix(network, comm_nodes=None, B=None):
    '''
    This function computes the modularity matrix
    for a specific group in the network.
    (a.k.a., generalized modularity matrix)

    Specifically,
    .. math::
        B^g_{i,j} = B_ij - \delta_{ij} \sum_(k \in g) B_ik
        m = \abs[\Big]{E}
        B_ij = A_ij - \dfrac{k_i k_j}{2m}
        OR...
        B_ij = \(A_ij - \frac{k_i^{in} k_j^{out}}{m}

    When `comm_nodes` is None or all nodes in `network`, this reduces to :math:`B`

    Parameters
    ----------
    network : nx.Graph or nx.DiGraph
        The network of interest
    comm_nodes : iterable (list, np.array, or tuple)
        List of nodes that defines a community
    B : np.matrix
        Modularity matrix of `network`

    Returns
    -------
    np.matrix
        The modularity of `comm_nodes` within `network`
    '''

    if comm_nodes is None:
        comm_nodes = list(network)
        return get_base_modularity_matrix(network)

    if B is None:
        B = get_base_modularity_matrix(network)

    # subset of mod matrix in g
    indices = [list(network).index(u) for u in comm_nodes]
    B_g = B[indices, :][:, indices]
    #print 'Type of `B_g`:', type(B_g)

    # B^g_(i,j) = B_ij - δ_ij * ∑_(k∈g) B_ik
    # i, j ∈ g
    B_hat_g = np.zeros((len(comm_nodes), len(comm_nodes)), dtype=float)

    # ∑_(k∈g) B_ik
    B_g_rowsum = np.asarray(B_g.sum(axis=1))[:, 0]
    if type(network) == nx.Graph:
        B_g_colsum = np.copy(B_g_rowsum)
    elif type(network) == nx.DiGraph:
        B_g_colsum = np.asarray(B_g.sum(axis=0))[0, :]

    for i in range(B_hat_g.shape[0]):
        for j in range(B_hat_g.shape[0]):
            if i == j:
                B_hat_g[i,j] = B_g[i,j] - 0.5 * (B_g_rowsum[i] + B_g_colsum[i])
            else:
                B_hat_g[i,j] = B_g[i,j]

    if type(network) == nx.DiGraph:
        B_hat_g = B_hat_g + B_hat_g.T

    return B_hat_g

def largest_eig(A):
    '''
        A wrapper over `scipy.linalg.eig` to produce
        largest eigval and eigvector for A when A.shape is small
    '''
    vals, vectors = eig(A)
    real_indices = [idx for idx, val in enumerate(vals) if not bool(val.imag)]
    vals = [vals[i].real for i in range(len(real_indices))]
    vectors = [vectors[i] for i in range(len(real_indices))]
    max_idx = np.argsort(vals)[-1]
    return np.asarray([vals[max_idx]]), np.asarray([vectors[max_idx]]).T

