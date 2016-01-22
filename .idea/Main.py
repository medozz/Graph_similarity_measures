"Author: Dezso Modos"
import igraph
import scipy
import numpy as np

"""
Description:
The following script will bring a graph similarity matrixs bsed on particular nodes and there selected neighbours.
Graphs should be in .ncol format or any other format, which igraph reads. Nodes in this case proteins, but any entities
can be.
Requieres a list of nodes, which have some kind of relatedness. Theese can be any kind of relatedness: compound biding,
cancer relatedness etc. This is a list of entities and a weight. The weight definitly is 1, if otherwiese not specified.
Output:
The script will give you a particular matrix of similarity between two graph, where the nodes are related.
"""



