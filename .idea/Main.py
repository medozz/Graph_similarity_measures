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

def import_nodes(file_name,separator,default_weight, header):
    """
    :param file_name: imported file name
    :param separator: separator of the import file
    :param default_weight: 0 if there is no default wheight and use the files wheights
    :param header: if true there is a header, if false there is no header
    :return: id_wieght dictionarry with string key and float weights
    """
    id_weight = {}
    inp = open(file_name)
    if header == True:
        header = inp.readline()
    for line in inp:
        line = line.split(sep)
        if len(line) > 1 and default_weight==0:
            id_weight[str(line[0])] = float(line[1])
        elif float(default_weight) > 0:
            id_weight[str(line[0])] = float(default_weight)
        else:
            id_weight[str(line[0])] = 1.0
    return id_weight


def giancomponenet(G):
    '''
    Simple function. Makes form a graph a giant componenet graph.
    '''
    print "Making giant componenet"
    components = G.clusters(IG.WEAK)
    gr2="giantcomponent.txt"
    igraph.save(components.giant(), gr2, format="ncol")
    G = igraph.Graph.Read_Ncol(open(gr2,"rb"),names=True, weights="if_present", directed=True)
    return G


def neighbors_flow_propagation(vertexes, information, used_vertexes, outgraf):
    """

    :param vertexes: The vertxes, whose neighbors are searched for
    :param information: The amount of informatioon which cvommes in
    :param used_vertexes: the set of used vertexies
    :return: added used vertexes, nighbor vertexes, for the next step
    """
    for vertex in vertexes:
        neighbors_of_vertex=vertex.neighbors()
        for one in neighbors_of_vertex:
            if one not in used_vertexes:
                one["Reach"]= one["Reach"]+information*0.5
    for one in vertexes.neighbors():
        used_vertexes.add(one)
    neighbor_vertexes = vertexes.neighbors()
    return used_vertexes, neighbor_vertexes, outgraf


def relatedness_count(ourgraf,id_weight,neighborhood_number):
    """
    :param ourgraf: igraph graph file
    :param id_weight: id_weights from the funcition import_nodes
    :param neighborhood_number int, tells how far4 should igraph chech the neighbors
    :return: igraph graph object with id_weight in

    Weigt propagation is not changable jet. Curently used the function as 2^-n*(income weight)
    Later can be used the edge weight of the particular graph, or thinking about different function.
    """
    related=[]
    FN=0
    ourgraf.vs["Reach"]=0
    for vertex in ourgraf.vs:
        if vertex['name'] in id_weight:
            vertex["Reach"] = vertex["Reach"]+id_weight[vertex['name']]
            used_vertexes=set(vertex)
            neighbor_vertexes=set(vertex)
            k=1
            while k<neighborhood_number+1:
                infromation = id_weight[vertex['name']]*2^-k
                used_vertexes, neighbor_vertexes, outgraf=neighbors_flow_propagation(neighbor_vertexes, infromation,
                                                                            used_vertexes, outgraf)
                k=k+1
    return ourgraf


def outwirte(outgraf, file_name, sep):
    out = open(file_name,"wb")
    out.write("UPID"+sep+"Reach"+"\n")
    for vertex in outgraf.vs:
        out.write(vertex['name']+"\t"+str(vertex["Reach"])+"\n")
    out.close()

#Running commands
id_weights = import_nodes(sample_nodes,"\t",1)
G = igraph.Graph.Read_Ncol(open("sample.csv","rb"),names=True, weights="if_present", directed=True)
G = giancomponenet(G)
G = relatedness_count(G, id_weights, 1)
outwirte(G, "sample_out.txt", "\t")
