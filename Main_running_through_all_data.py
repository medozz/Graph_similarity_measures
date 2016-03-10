import Main as M
import Cell_line_eater
import os
import time
import read_in_dan_descriptor_make_sam_from_graphs as Dan_fromat
#Globals:
graph = "ominpath.ncol"#used graph
folder = "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/" #used folder
sample_file = "dan_combination_train_signalink_targets_and_neighbours.csv" # To make accoruate file format using this file
Krishna_file = "cell_line_tissue_dictionary.txt" #Tissue translating if something is missing than this is
#  corrigated based here the tissue annotation. Curently after the graph based test.

#Import annotations
SwissProtUniProtIds = M.uniprotin(folder+"uniprot-homo+sapiens.tab")
gene_name_uniprot_library = M.chip_annotation_1_to_chip_annotation_2(folder+"GPL13667-15572_annotation.csv", 14, 21, SwissProtUniProtIds)
print "Annotations are ready they are in the memory"

#First line creats the expressions of the cell lines. Theese are the so called .expr files.
# Theese are the list of expressing mRNAs in the cell lines. The graph vbased analyisis is always this is the first
# step. The function have two type of argument  1:
# SD as -1 SD expression do not counted as expressed
# PT as Percentage of expression less than 25 counted as not expressed.
# This will be first SD based and than percentage based.
Cell_line_eater.cell_line_extractor(folder, "gex.csv", "SD")
print "Cell lines data are created. "
expressions=[]
for each in os.listdir(folder):
    if each.endswith("SD.expr"):
        expressions += [each]

def run_through_expression(neigbourhood, propagation_type, expressions):
    """
    This function creates the cell line graph based descriptions according to the paramters
    :param neigbourhood: int, tells how far should igraph chek the neighbours
    :param propagation_type: 1 if a vertex is reached from an oter one added to used vertexes
    2 add to used vertexes only the end of the cycle to prevent infomration backflow.
    3 not add to used vertexes, allow backword infomration flow
    :param expressions: expression files containing list
    :return:
    """
    for each in expressions:
        M.graph_from_expression_file_graph(folder+each, graph, gene_name_uniprot_library, neigbourhood, propagation_type)
        b = time.clock()
        print "Cell line completed:", each, "Time ellapsed since start:", (b-a)/60, "minutes"
        print each, "done"

print "SD Expressions are Done :)"
