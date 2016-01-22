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
    :return: id_wieght list of tuples with string id and float weight pairs
    """
    id_weight = []
    inp = open(file_name)
    if header == True:
        header = inp.readline()
    for line in inp:
        line = line.split(sep)
        if len(line) > 1 and default_weight==0:
            id_weight.append((str(line[0]),float(line[1])))
        elif float(default_weight) > 0:
            id_weight.append((str(line[0]),float(line[1])))
        else:
            id_weight.append((str(line[0]),1.0))
    return id_weight









def giancomponenet(G):
    '''
    Simple function. Makes form a graph a giant componenet graph.
    '''
    print "Making giant componenet"
    import igraph as IG
    components = G.clusters(IG.WEAK)
    gr2="giantcomponent.txt"
    IG.save(components.giant(), gr2, format="ncol")
    G = IG.Graph.Read_Ncol(open(gr2,"rb"),names=True, weights="if_present", directed=True)
    return G


def CRfilereader(grafunk,inp):
    '''
    Reads cancer related data file.
    '''
    cancerrelated=[]
    inp=open(inp,"rb")
    FN=0
    grafunk.vs["CRtype"]=0
    for line in inp:
        line=line.strip()
        line=line.split("\t")
        try:
            int (line[1])
            if int(line[1])==1:
                cancerrelated.append(line[0])
        except:
            continue

    for a in grafunk.vs:

        if a['name'] in cancerrelated:
            a["CRtype"]=1
            szomszed=a.neighbors()
            for elem in szomszed:
                if elem["CRtype"]!=1:
                    elem["CRtype"]=2
                    FN=FN+1
    print "First neighbors: ", FN
    print "cancer related: ",len(cancerrelated)
    return grafunk


