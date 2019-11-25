#!/usr/bin/env python

# Explore the network

import networkx as nx
import matplotlib.pyplot as plt
from utils import get_modularity
import igraph as ig
import sys

#---------------------------------------------------------------------------
def hard_cluster(fname):
    G = ig.read(fname, format="gml")
    c = G.community_walktrap().as_clustering()
    return c
#---------------------------------------------------------------------------

if __name__ == "__main__":
    #fname = '../data/RectumMicrobiomeNetwork/stool_network.gml.network'
    pre_fname = '../data/RectumMicrobiome_PrePost/stool_network_pre.gml.txt.network'
    post_fname = '../data/RectumMicrobiome_PrePost/stool_network_post.gml.txt.network'

    #read the pre/post network
    pre_G = ig.read(pre_fname, format="gml")
    post_G = ig.read(post_fname, format="gml")

    #hard clusterig
    pre_c = pre_G.community_walktrap().as_clustering()
    post_c = post_G.community_walktrap().as_clustering()

    #print the informaiton of each cluster
    print(pre_c)
    print(post_c)

    #print the size of clusters
    print("pre network:")
    for cluster in pre_c:
        print(len(cluster))
    print("------------------------------------------------------------------")
    print("post network:")
    for cluster in post_c:
        print(len(cluster))

    #Modularity
    print("The modularity of the pre network: {}".format(pre_G.modularity(pre_c)))
    print("The modularity of the post network: {}".format(post_G.modularity(post_c)))




#--------------------------------------------------------------------------------------------
# else:
#     print("igraph not working")
#     # networkx
#     G = nx.read_gml(fname, label='id')
#
#     print(nx.info(G))
#
#     # Min weight: 0.670000184425068
#     # print(min([x['weight'] for _,_,x in G.edges(data=True)]))
#
#     from networkx.algorithms.community import greedy_modularity_communities
#
#     # Clauset-Newman-Moore greedy modularity maximization, modularity: 0.44539
#     c = greedy_modularity_communities(G)
#
#     community_assignment = {}
#
#     for i, community in enumerate(c):
#         for node in community:
#             community_assignment[node] = i
#
#     print(get_modularity(G, community_assignment))
