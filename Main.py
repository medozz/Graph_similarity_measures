"Author: Dezso Modos"
import igraph
import scipy
import numpy as np
import time
import os
import re
import string

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
The script only uses the giant componenet of the graph.
"""

def nodest_from_exp_file(exp_file, gene_name_to_uniprot):
    inp = open(exp_file)
    gene_name_nodeset = set()
    for line in inp:
        gene_name_nodeset.add(line.strip())
    inp.close()
    up_node_set = set()
    for node in gene_name_nodeset:
        for upid in gene_name_to_uniprot[node]:
            up_node_set.add(upid)
    return up_node_set


def prepare_graph(nodeset, graph_ncol, cell_line_graph):
    """
    This function reads the graph ncol file and search wheteher the particular nodes expressions are presented, or not.
    :param nodeset: the nodes ids in appropariate format (a set of ids)
    :param graph_ncol: the import ncol graph
    :param cell_line_graph: The return file name
    :return: it returns the porticular graph with the header of the inport file, but in ncol.
    """
    inp = open(graph_ncol)
    outedges=set()
    for edge in graph_ncol:
        edge = edge.split(" ")
        if edge[0].strip() and edge[1].strip() in nodeset:
            outedges.add(" ".join(edge).strip())
    inp.close()

    out = open(cell_line_graph, "wb")
    for edge in outedges:
        out.write(edge)+"\n"
    out.close()


def create_node_weight_file_from_gen_descritor(gene_name_uniprot_library, descriptorfile, separator, cell_line_column,
                                               gene_column_start, descriptortype, folder):
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
    log=open(folder+descriptortype+"_cellist_making_log.log", "wb")
    while k<len(header):
        if k>gene_column_start:
            gene_column_dico[header[k]] = k
        k = k+1
    failures=set()
    for line in inp:
        line = line.split(separator)
        out = open(folder+line[cell_line_column]+descriptortype+".celist","wb")
        for gene in gene_column_dico:
            try:
                for upid in gene_name_uniprot_library[gene]:
                    out.write(upid)
                    out.write("\t")
                    out.write(str(line[gene_column_dico[gene]])+"\n")
            except:
                failures.add(gene)
        out.close()

    log.write("\n"+str(len(failures))+"\n"+" ".join(failures))
    log.close()


def import_nodes(file_name, sep, default_weight, header):
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
            id_weight[str(line[0])] = abs(float(line[1]))
        elif float(default_weight) > 0:
            id_weight[str(line[0])] = float(default_weight)
        else:
            id_weight[str(line[0])] = 1.0
    return id_weight


def giancomponenet(G):
    """
    Simple function. Makes form a graph a giant componenet graph.
    """
    print "Making giant component"
    components = G.clusters(igraph.WEAK)
    gr2="giantcomponent.txt"
    igraph.save(components.giant(), gr2, format="ncol")
    G = igraph.Graph.Read_Ncol(open(gr2, "rb"), names=True, weights="if_present", directed=True)
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
        neighbors_of_vertex = vertex.neighbors()
        for one in neighbors_of_vertex:
            if one not in used_vertexes:
                one["Reach"] = one["Reach"]+information*0.5
            if one not in neighbor_vertexes:
                neighbor_vertexes.append(one)
            if type_ == 1:
                used_vertexes.append(one)
        if type_ == 2:
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
    Weights should be positive numbers. Negative numbers are not meaningful in this context. To use
    negative numbers the graph should be directed.
    """
    #a=float(time.clock())
    #FN=0
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
                infromation = abs(id_weight[vertex['name']])*(2**(-k))
                used_vertexes, neighbor_vertexes=neighbors_flow_propagation(neighbor_vertexes, infromation,
                                                                            used_vertexes, propagation_type)
                k=k+1
            """
            if counter%100==0:
                print "From", all_sources, "in the graph",100*counter/all_sources, "%"
                b=float(time.clock())
                print b-a
            """
    print "The primarry uniprot IDS in the graph:", counter, "Percentage:", (float(counter)/float(all_sources))*100
    return ourgraf


def outwirte(outgraf, file_name, sep):
    out = open(file_name,"wb")
    out.write("UPID"+sep+"Reach"+"\n")
    for vertex in outgraf.vs:
        out.write(vertex['name']+"\t"+str(vertex["Reach"])+"\n")
    out.close()

def chip_annotation_1_to_chip_annotation_2(chip_annotation_file, id1_col, id2_col, swissprott_uniprot_ids):
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
                        if id2 in swissprott_uniprot_ids:
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


def uniprotin(uniprotfile):
    inp = open(uniprotfile)
    inp.readline()
    SwissProtUniProtIds=set()
    for line in inp:
        line=line.split("\t")
        SwissProtUniProtIds.add(line[0])
    return SwissProtUniProtIds


#Running commands

folder = "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/run/"
graph = folder+"Signor_2016_02_05.ncol"

SwissProtUniProtIds=uniprotin(folder+"uniprot-homo+sapiens.tab")

gene_name_uniprot_library = chip_annotation_1_to_chip_annotation_2(folder+"GPL13667-15572_annotation.csv", 14, 21, SwissProtUniProtIds)
create_node_weight_file_from_gen_descritor(gene_name_uniprot_library,
                                           folder+"cell_line_gene_distance_fingerprints.csv", ",",
                                           1, 2, "cell_line_gene_distance_affy_translation_only_SP", folder)


results = []
for each in  os.listdir(folder):
    if each.endswith(".celist"):
         results += [each]
print "results:", results

expressions=[]
for each in  os.listdir(folder):
    if each.endswith("SD.expr"):
        expressions += [each]
"""
trues=[]
not_trues=[]
for each in expressions:
    each2=each.replace("_SD.expr", "cell_line_gene_distance_affy_translation_only_SP.celist")
    if each2 in results:
        trues.append(each2)
    else:
        not_trues.append((each, each2))
"""
def graph_from_expression_file_graph(expression_file, graph_file, GENE_name_uniprot):
    expression_set=nodest_from_exp_file(expression_file, GENE_name_uniprot)

    prepare_graph(expression_set,graph_file,expression_file.reaplace("expr", "ncol"))
    id_weights = import_nodes(folder+cell_line,"\t",0, 0)
    expression_graph_file = open(expression_file.reaplace("expr", "ncol"))
    G = igraph.Graph.Read_Ncol(expression_graph_file,names=True, weights="if_present", directed=True)
    G = giancomponenet(G)
    G = relatedness_count(G, id_weights, 2, 1) #according to Krishna Neighborhood will be 2 propagation type will be 1
    cell_line = expression_file.replace("_SD.expr", "cell_line_gene_distance_affy_translation_only_SP.celist")
    outwirte(G, string.replace(folder+cell_line, ".celist", "Signor_no_backward_propagation_three_neighbor.celdesc"), "\t")
    # Line above should be rewritten at any paramters run

a=float(time.clock())
for each in expressions:
    graph_from_expression_file_graph(folder+each, graph, gene_name_uniprot_library)
    b=time.clock()
    print "Cell line completed:", each, "Time ellapsed since start:", (b-a)/60, "minutes"
    print each, "done"
print "Done :)"
gene_name_uniprot_library=""
"""
print len(expressions), len(results), len (trues)
for a in not_trues:
    print a

#print gene_name_uniprot_library

a=float(time.clock())
#graph = open("Reactome_2016_01_22_onlySP.ncol")
graph = open("UP_string_2016_01_22.ncol")
id_weights = import_nodes("NCI-SNU-16cell_line_strength_new_translation.celist","\t",0, 0)
G = igraph.Graph.Read_Ncol(open("UP_string_2016_01_22.ncol","rb"),names=True, weights="if_present", directed=True)
G = giancomponenet(G)
G = relatedness_count(G, id_weights, 2, 1) #according to Krishna Neighborhood will be 2 propagation type will be 1
outwirte(G, "STRING_NCI-SNU-16cell_line.celdesc", "\t")
b=time.clock()
print " New upit translation, Time ellpased since start:", b-a

a=float(time.clock())
#graph = open("Reactome_2016_01_22_onlySP.ncol")
graph = open("UP_string_2016_01_22.ncol")
id_weights = import_nodes("NCI-SNU-16cell_line_strength.celist","\t",0, 0)
G = igraph.Graph.Read_Ncol(open("UP_string_2016_01_22.ncol","rb"),names=True, weights="if_present", directed=True)
G = giancomponenet(G)
G = relatedness_count(G, id_weights, 2, 1) #according to Krishna Neighborhood will be 2 propagation type will be 1
outwirte(G, "STRING_NCI-SNU-16cell_line.celdesc", "\t")
b=time.clock()
print "Time ellpased since start:", b-a
"""