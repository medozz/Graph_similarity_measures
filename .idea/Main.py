"Author: Dezso Modos"
import igraph
import scipy
import numpy as np
import time
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
    inp = open(descriptorfile)
    header = inp.readline()
    header = header.split(separator)
    k=0
    gene_column_dico={}
    while k<len(header):
        if k>gene_column_start:
            gene_column_dico[header[k]] = k
        k = k+1
    failures=set()
    for line in inp:
        line = line.split(separator)
        out = open("/home/dm729/PycharmProjects/Graph_similarity_measures/"+line[cell_line_column]+descriptortype+".celist","wb")
        for gene in gene_column_dico:
            try:
                for upid in gene_name_uniprot_library[gene]:
                    out.write(upid)
                    out.write("\t")
                    out.write(str(line[gene_column_dico[gene]])+"\n")
            except:
                failures.add(gene)

        out.close()
    print len(failures), failures


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
    a=float(time.clock())
    related=[]
    FN=0
    counter=0
    all_sources = len(id_weight)
    ourgraf.vs["Reach"]=0
    for vertex in ourgraf.vs:
        if vertex['name'] in id_weight:
            vertex["Reach"] = vertex["Reach"]+id_weight[vertex['name']]
            used_vertexes=[vertex]
            neighbor_vertexes=[vertex]
            k=0
            counter = counter+1
            while k<neighborhood_number:
                infromation = id_weight[vertex['name']]*(2**(-k))
                used_vertexes, neighbor_vertexes=neighbors_flow_propagation(neighbor_vertexes, infromation,
                                                                            used_vertexes, propagation_type)
                k=k+1
            if counter%100==0:
                print "From", all_sources, 100*counter/all_sources, "%"
                b=float(time.clock())
                print b-a
    return ourgraf


def outwirte(outgraf, file_name, sep):
    out = open(file_name,"wb")
    out.write("UPID"+sep+"Reach"+"\n")
    for vertex in outgraf.vs:
        out.write(vertex['name']+"\t"+str(vertex["Reach"])+"\n")
    out.close()

def chip_annotation_1_to_chip_annotation_2(chip_annotation_file, id1_col, id2_col):
    """
    Reads in an Affymetrix chip annotation file from GEO and gives a dictionarry, which contains id1: set(id2s)
    The affymetrix files almost any properties are hardcoded, like the separators between and within columns.
    :param chip_annotation_file: Affymetrix annoatation file of the chip
    :param id1_col: source column
    :param id2_col: target identifier column
    :return: id1 to set(id2s) dictionarry
    """
    inp = open(chip_annotation_file)
    gene_name_to_uniprot_dictionarry={}
    for line in inp:
        if line[0] != "#":
            line = line.split("\t")
            for id1 in line[id1_col].split(" /// "):
                if id1 in gene_name_to_uniprot_dictionarry:
                     for id2 in line[id2_col].split(" /// "):
                         gene_name_to_uniprot_dictionarry[id1].add(id2)
                else:
                    gene_name_to_uniprot_dictionarry[id1]=set()
                    for id2 in line[id2_col].split(" /// "):
                         gene_name_to_uniprot_dictionarry[id1].add(id2)
    return gene_name_to_uniprot_dictionarry


def read_uniprot_dictionarry(up_file,sep,reviewed_col):
    """
    REads uniprot file to make gene name to uniprot ID
    :param up_file: uniprot file
    :param sep: sepoarator as in the downloaded file
    :param reviewed_col: place of the reviewed col
    :return: gene_name_up_dic: gene_name set(upid)
    """
    inp=open(up_file)
    gene_name_up_dic = {}
    for line in inp:
        line = line.split(sep)
        if line[reviewed_col]== "reviewed":
            if line[0] not in gene_name_up_dic:
                gene_name_up_dic[line[0]] = set()
                gene_name_up_dic[line[0]].add(line[1])
            else:
                gene_name_up_dic[line[0]].add(line[1])
    return gene_name_up_dic

#Running commands

gene_name_uniprot_library = chip_annotation_1_to_chip_annotation_2("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/GPL13667-15572_annotation.csv", "\t", 3)
create_node_weight_file_from_gen_descritor(gene_name_uniprot_library,
                                           "/home/dm729/PycharmProjects/Graph_similarity_measures/cell_line_gene_distance_fingerprints.csv", ",",
                                           1,2,"cell_line_strength_new_translation")

a=float(time.clock())
graph = open("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/Reactome_2016_01_22_onlySP.ncol")
id_weights = import_nodes("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/NCI-SNU-16cell_line_strength_new_translation.celist","\t",0, 0)
G = igraph.Graph.Read_Ncol(open("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/Reactome_2016_01_22_onlySP.ncol","rb"),names=True, weights="if_present", directed=True)
G = giancomponenet(G)
G = relatedness_count(G, id_weights, 2, 1) #according to Krishna Neighborhood will be 2 propagation type will be 1
outwirte(G, "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/Reactome_NCI-SNU-16cell_line.celdesc", "\t")
b=time.clock()
print " New upit translation, Time ellpased since start:", b-a

a=float(time.clock())
graph = open("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/Reactome_2016_01_22_onlySP.ncol")
id_weights = import_nodes("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/NCI-SNU-16cell_line_strength.celist","\t",0, 0)
G = igraph.Graph.Read_Ncol(open("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/Reactome_2016_01_22_onlySP.ncol","rb"),names=True, weights="if_present", directed=True)
G = giancomponenet(G)
G = relatedness_count(G, id_weights, 2, 1) #according to Krishna Neighborhood will be 2 propagation type will be 1
outwirte(G, "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/Reactome_NCI-SNU-16cell_line.celdesc", "\t")
b=time.clock()
print "Time ellpased since start:", b-a