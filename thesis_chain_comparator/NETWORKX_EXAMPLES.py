"""
Here I demonstrate how to assess the relationship between networks using
the NetworkX library. More specifically, I am comparing whether nodes
in one network for one edge type are shared in another network whose
edges are of a different identity.
"""

import networkx as nx
from numpy import diag, linalg

# a call to nx.draw(G) can be made to view the graphs

# --------------------------------------------------------------------------------
# // no relationship //

# the base graph
G1 = nx.Graph()
G2 = nx.Graph()

# add nodes and edges
G1.add_nodes_from([1, 2, 3, 4, 5, 6])
G1.add_edges_from([(1, 2), (2, 3), (3, 4)])

G2.add_nodes_from([1, 2, 3, 4, 5, 6])
G2.add_edges_from([(5, 6)])

# the adjacencies
A1 = nx.adj_matrix(G1).todense()
A2 = nx.adj_matrix(G2).todense()

D = diag(A1 * A2)
print('No relationship: ', D)

# --------------------------------------------------------------------------------
# // direct superimposition //

# the base graph
G1 = nx.Graph()
G2 = nx.Graph()

# add nodes and edges
G1.add_nodes_from([1, 2, 3, 4])
G1.add_edges_from([(1, 2), (2, 3), (3, 4)])

G2.add_nodes_from([1, 2, 3, 4])
G2.add_edges_from([(2, 3)])

# the adjacencies
A1 = nx.adj_matrix(G1).todense()
A2 = nx.adj_matrix(G2).todense()

# ** rank-nullity theorem example
G1G2 = A1 * A2  # match thesis
DIM = G1G2.shape[1]
RANK = linalg.matrix_rank(G1G2)
NULLITY = DIM - RANK  

D = diag(G1G2)
print('Direct superimposition: ', D, 'Nullity: ', NULLITY)

# --------------------------------------------------------------------------------
# // indirect superimposition //

# the base graph
G1 = nx.Graph()
G2 = nx.Graph()

# add nodes and edges
G1.add_nodes_from([1, 2, 3, 4, 5])
G1.add_edges_from([(1, 2), (2, 3), (3, 4)])

G2.add_nodes_from([1, 2, 3, 4, 5])
G2.add_edges_from([(1, 5)])

# the adjacencies
A1 = nx.adj_matrix(G1).todense()
A2 = nx.adj_matrix(G2).todense()

D = diag(A1 * A2)
print('Indirect superimposition: ', D)

# --------------------------------------------------------------------------------
# // pseudomembership //

# the base graph
G1 = nx.Graph()
G2 = nx.Graph()

# add nodes and edges
G1.add_nodes_from([1, 2, 3, 4])
G1.add_edges_from([(1, 2), (3, 4)])

G2.add_nodes_from([1, 2, 3, 4])
G2.add_edges_from([(1, 2), (3, 4)])

# the adjacencies
A1 = nx.adj_matrix(G1).todense()
A2 = nx.adj_matrix(G2).todense()

# ** rank-nullity theorem example
G1G2 = A1 * A2  # match thesis
DIM = G1G2.shape[1]
RANK = linalg.matrix_rank(G1G2)
NULLITY = DIM - RANK

D = diag(A1 * A2)
print('Pseudomembership: ', D, 'Nullity: ', NULLITY)

