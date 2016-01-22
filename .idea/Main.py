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
According to Krishnaq the inputs will be weighted, but kept as simple as possible, so no back propagation and forward
propagation is used. For weights the different expression changes, methylations are used. For mutations  only the
dreiver mutations used.
"""
def create_node_weight_file_from_gen_descritor(gene_name_uniprot_library, descriptorfile, separator, cell_line_column,
                                               gene_column_start, descriptortype):
    """
    This function reads in a gene based descritor and forms a specific format for further annotation
    :param gene_name_uniprot_library: dictionarry for gene name to set(uniprot_id) - I can not be sure about this
    :param descriptorfile: descriptorfile made by Krishna
    :param separator: the separator between the columns
    :param cell_line_column: clumn of cell line
    :return: a cell line named descritor with values
    """
    inp=oepn(descriptorfile)
    header = inp.readline()
    header = header.split(separator)
    k=0
    gene_column_dico={}
    while k<len(header):
        if k>gene_column_start:
            gene_column_dico[header[k]] = k
        k = k+1
    for line in inp:
        line = line.split(separator)
        out = open(line[cell_line_column]+descriptortype+".celist")
        for gene in gene_column_dico:
            for upid in gene_name_uniprot_library[gene]:
                out.write(upid)
                out.write("\t")
                out.write(str(line[gene_column_dico[gene]]))
        out.close()


def import_nodes(file_name,sep,default_weight, header):
    """
    :param file_name: imported file name
    :param sep: separator of the import file
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
    components = G.clusters(igraph.WEAK)
    gr2="giantcomponent.txt"
    igraph.save(components.giant(), gr2, format="ncol")
    G = igraph.Graph.Read_Ncol(open(gr2,"rb"),names=True, weights="if_present", directed=True)
    return G


def neighbors_flow_propagation(vertexes, information, used_vertexes, type_):
    """

    :param vertexes: The vertxes, whose neighbors are searched for
    :param information: The amount of informatioon which cvommes in
    :param used_vertexes: the set of used vertexies
    :param type_: 1 if a vertex is reached from an oter one added to used vertexes
    2 add to used vertexes only the end of the cycle to prevent infomration backflow.
    3 not add to used vertexes, allow backword infomration flow
    :return: added used vertexes, nighbor vertexes, for the next step
    """
    neighbor_vertexes = []
    if type_==3:
        used_vertexes = []
    for vertex in vertexes:
        neighbors_of_vertex=vertex.neighbors()
        for one in neighbors_of_vertex:
            if one not in used_vertexes:
                one["Reach"]= one["Reach"]+information*0.5
            if one not in neighbor_vertexes:
                neighbor_vertexes.append(one)
            if type_==1:
                used_vertexes.append(one)
        if type_==2:
            for vertex in neighbors_of_vertex:
                used_vertexes.append(vertex)
    return used_vertexes, neighbor_vertexes


def relatedness_count(ourgraf,id_weight,neighborhood_number, propagation_type):
    """
    :param ourgraf: igraph graph file
    :param id_weight: id_weights from the funcition import_nodes
    :param neighborhood_number int, tells how far4 should igraph chech the neighbors
    :param propagation_type 1 if a vertex is reached from an oter one added to used vertexes
    2 add to used vertexes only the end of the cycle to prevent infomration backflow.
    3 not add to used vertexes, allow backword infomration flow
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
            used_vertexes=[vertex]
            neighbor_vertexes=[vertex]
            k=0
            while k<neighborhood_number:
                infromation = id_weight[vertex['name']]*(2**(-k))
                used_vertexes, neighbor_vertexes=neighbors_flow_propagation(neighbor_vertexes, infromation,
                                                                            used_vertexes, propagation_type)
                k=k+1
    return ourgraf


def outwirte(outgraf, file_name, sep):
    out = open(file_name,"wb")
    out.write("UPID"+sep+"Reach"+"\n")
    for vertex in outgraf.vs:
        out.write(vertex['name']+"\t"+str(vertex["Reach"])+"\n")
    out.close()

#Running commands
a = open("/home/dm729/PycharmProjects/Graph_similarity_measures/sample.csv")
id_weights = import_nodes("/home/dm729/PycharmProjects/Graph_similarity_measures/sample_weights","\t",0, 0)
G = igraph.Graph.Read_Ncol(open("/home/dm729/PycharmProjects/Graph_similarity_measures/sample.csv","rb"),names=True, weights="if_present", directed=True)
G = giancomponenet(G)
G = relatedness_count(G, id_weights, 2, 3)
outwirte(G, "/home/dm729/PycharmProjects/Graph_similarity_measures/sample_out.txt", "\t")
