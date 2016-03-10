import os
import numpy as np
import string

def read_in_dan_descriptor(dan_descriptor_file):
    """
    Reads in a descriptor file from Dan and returns a list of th descriptor names in a tupple.
    The first element of the tupple is a
    :param dan_descriptor_file: Dan descriptor file in accepted format.
    :return: list of tuples as described above.
    """
    inp=open(dan_descriptor_file)
    returnlist=[]
    inp.readline() #header
    for line in inp:
        line = line.strip()
        line = line.split(",")
        cell_line = line[0].split(">")[1]
        returnlist.append((line[0], cell_line))
    inp.close()
    return returnlist

def import_data_from_cell_line_file(cell_line_file, sep):
    """
    This function will import the descriptor from the cell line file to add a dictionarry.
    Important it also normalieze the values to have a better value in the model.
    :param cell_line_file: Imported cell line descriptor
    :param sep: The separator between gene name and value
    :return: the dictionarry containing Gene:value pairs
    """
    inp = open(cell_line_file)
    inp.readline()

    list_of_values = []
    for line in inp:
        line = line.strip()
        line = line.split(sep)
        list_of_values.append(float(line[1]))
    inp.close()
    mean = np.mean(list_of_values)
    inp = open(cell_line_file)
    inp.readline()
    return_dico={}
    for line in inp:
        line = line.strip()
        line = line.split(sep)
        return_dico[line[0]] = np.float32(float(line[1])/mean)
    return return_dico

def construct_descriptors_for_additional_cell_lines(krishna_translation_file, krishna_file_separator, folder, descriptor_description):
    """
    This function will construct of the not presented cell lines (or whathever)descriptors to the average of a
    descriptor category. In this case this is the cell lines which are from the same tissue.
    :param krishna_translation_file: this file contains some kind of dictionarry cell line: tissue (smaller to larger
     category)
    :param folder: the folder where are the descriptors
    :param descriptor_description: the descriptors file external. Here is usally the methods, the graph or any paramter.
    Please be aveare to use all the file extension after the cell line (small category), or the script will not work.
    :return: In the directory there will be average results of the lacking cell lines.
    """
    dictionarry_tissue_to_cell_line_set = {}
    dictionarry_cell_line_to_tissue = {}
    sep = krishna_file_separator
    inp = open(krishna_translation_file)
    cell_lines_all_in_dream = set()
    inp.readline() #Header
    for line in inp:
        line = line.strip()
        line = line.split(sep)
        if line[0] not in dictionarry_cell_line_to_tissue:
            dictionarry_cell_line_to_tissue[line[0]] = line[1]
        else:
            print line, "Failure check Krishan's source file!"
        if line[1] not in dictionarry_tissue_to_cell_line_set:
            dictionarry_tissue_to_cell_line_set[line[1]]=set()
        dictionarry_tissue_to_cell_line_set[line[1]].add(line[0])
        cell_lines_all_in_dream.add(line[0])

    print "All cell lines in Dream:",  cell_lines_all_in_dream

    cell_descriptor_files = []
    for each in os.listdir(folder):
        if each.endswith(descriptor_description):
            cell_descriptor_files.append(each)

    cell_lines_with_descriptors = set()
    for file_name in cell_descriptor_files:
        cell_line = string.replace(file_name, descriptor_description, "")
        cell_lines_with_descriptors.add(cell_line)

    all_cell_lines_dictionarry = {}
    genelist=[]
    for cell_line in cell_lines_with_descriptors:
        cell_line_file = folder+cell_line+descriptor_description
        cell_line_dictionarry = import_data_from_cell_line_file(cell_line_file, "\t")
        all_cell_lines_dictionarry[cell_line] = cell_line_dictionarry

        for gene in cell_line_dictionarry:
            if gene not in genelist:
                genelist.append(gene)

    no_descriptor_cell_lines = cell_lines_all_in_dream - cell_lines_with_descriptors
    out_no_information = open(folder+"Noinfroamtion_in_"+string.replace(descriptor_description,".celdesc","log"), "wb")

    print"No descriptor cell lines:", no_descriptor_cell_lines

    for no_information_cell_line in no_descriptor_cell_lines:
        print no_information_cell_line
        out_no_information.write(no_information_cell_line+"\n")
        gene_value_pair = {}
        for gene in genelist:
            gene_value_pair[gene]={}
            gene_value_pair[gene]["values"] = []

        for cell_line in dictionarry_tissue_to_cell_line_set[dictionarry_cell_line_to_tissue[no_information_cell_line]]:
            if cell_line != no_information_cell_line:
                for gene in genelist:
                    if gene in all_cell_lines_dictionarry[cell_line]:
                        gene_value_pair[gene]["values"].append(float(all_cell_lines_dictionarry[cell_line][gene]))
                    else:
                        gene_value_pair[gene]["values"].append(0)

        out = open(folder+no_information_cell_line+descriptor_description, "wb")
        out.write("UPID\tReach")
        for gene in genelist:
            gene_value_pair[gene]["mean"] = np.mean(gene_value_pair[gene]["values"])
            gene_value_pair[gene]["median"] = np.median(gene_value_pair[gene]["values"])
            gene_value_pair[gene]["SD"] = np.std(gene_value_pair[gene]["values"])
            out.write(gene+"\t"+str(gene_value_pair[gene]["mean"])+"\n")
        out.close()
    out_no_information.close()


def open_cell_lines_write_out_descriptor(folder, Dan_format_list_of_cell_lines, descriptor_description, outputfile):
    """
    :param folder: Place of cell descriptors
    :param Dan_format_list_of_cell_lines: List of cell descriptors in Dan file
    :param descriptor_description: The external descripion of the descriptor, which is not the cell line.
    :param outputfile: The return file_the Graph descriptor
    :return:
    """
    cell_lines=set()
    for things in Dan_format_list_of_cell_lines:
        #if things[1] != "MDA-MB-175-VII" and things[1] != "NCI-H1437":
        cell_lines.add(things[1])
    all_cell_lines_dictionarry = {}
    genelist=[]
    for cell_line in cell_lines:
        cell_line_file = folder+cell_line+descriptor_description
        cell_line_dictionarry = import_data_from_cell_line_file(cell_line_file, "\t")
        all_cell_lines_dictionarry[cell_line] = cell_line_dictionarry

        for gene in cell_line_dictionarry:
            if gene not in genelist:
                genelist.append(gene)
    out=open(outputfile, "wb")
    out.write("ID")
    for gene in genelist:
        out.write("," + gene)
    out.write("\n")
    for cell_line_compound_pair in Dan_format_list_of_cell_lines:
        out.write(cell_line_compound_pair[0])

        for gene in genelist:
            #if cell_line_compound_pair[1] != "MDA-MB-175-VII" \
            #        and cell_line_compound_pair[1] != "NCI-H1437":
            if gene in all_cell_lines_dictionarry[cell_line_compound_pair[1]]:
                out.write(","+str(all_cell_lines_dictionarry[cell_line_compound_pair[1]][gene]))
            else:
                out.write(",0")
        out.write("\n")
    out.close()

sample_file="dan_combination_train_signalink_targets_and_neighbours.csv"
Krishna_file="cell_line_tissue_dictionary.txt"
folder ="/scratch/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/run/"
graph_descritor_type = "cell_line_gene_distance_affy_translation_only_SPSignor_no_backward_propagation_DE_second_neighbor.celdesc"
gene_expression_descriptor = "cell_line_gene_distance_affy_translation_only_SP.celist"
final_file_name = "cell_differetnial_expression.csv"#dezso_combination_test_Signor_2nd_neighbors_GeneExpression_Mean_minus_one_SD_not_exp_mean"
construct_descriptors_for_additional_cell_lines(folder+Krishna_file, "\t", folder, gene_expression_descriptor)
dan_descriptor_file = read_in_dan_descriptor(folder+sample_file)
open_cell_lines_write_out_descriptor(folder, dan_descriptor_file, gene_expression_descriptor,
                                     final_file_name)


