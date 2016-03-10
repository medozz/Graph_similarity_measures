"Author: Dezso Modos"
"""
This script will eat the cell lines expressional data and calculate the not expressed genes. Not expressioncan be
calculated by two way:
1. RMA mean expression in the cell line and than -1 SD is calculated not expressed.
2. Below 25 percentile of expression calculated not expressed.
The script returns a "cellline_name.expr" file, which will be eaten by main. It is a file containing the gene names
of each expressed genes. Main will be rewritten in graph analysis to incorporate only theese genes.
"""
"!TO_DO: Rewrite main"

import numpy as np
import string


def calculate_expressed_percentage_based(values, percentile):
    """
    Returns from the values list the given percentile
    """
    return float(np.percentile(values, float(percentile)/100.0))


def calculate_expressed_sd_based(values, SDs):
    """
    Calculate the mean-SDs*SD value and returns it.
    """
    mean_ = np.mean (values)
    SD_ = np.std(values)
    return float(mean_)-float(SDs)*float(SD_)


def cell_line_extractor(folder ,cell_line_expression_file, method):
    """
    Thias is the main function for this work. The separators are coded.
    :param cell_line_expression_file: Dream challange particular cell line file
    :param method: The method falg. SD it calculates less than 1 SD gewnes as not expressed. PT it calculates less than
    25 percentage ase not expressed
    :return: Little files containing the names of genes which are presented in a particular cell line. The extension is
    .expr
    """
    inp = open(folder+cell_line_expression_file)
    header = inp.readline()
    cell_lines = header.split(",")[1:]
    cell_line_dictionarry = {}
    for cell_line in cell_lines:
        cell_line = (cell_line.replace('"', '')).strip()
        cell_line_dictionarry[cell_line] = {}
        cell_line_dictionarry[cell_line]["values"] = []
    for line in inp:
        line = line.split(",")
        k = 0
        while k< len(cell_lines):
            cell_line = (cell_lines[k].replace('"', '')).strip()
            cell_line_dictionarry[cell_line][line[0]] = float(line[k+1])
            cell_line_dictionarry[cell_line]["values"].append(float(line[k+1]))
            k = k+1

    for cell_line in cell_line_dictionarry:
        if method == "SD":
            critical_value = calculate_expressed_sd_based(cell_line_dictionarry[cell_line]["values"], 1)
            cell_line_dictionarry[cell_line]["critical_value"] = critical_value
        if method == "PT":
            critical_value = calculate_expressed_percentage_based(cell_line_dictionarry[cell_line]["values"], 25)
            cell_line_dictionarry[cell_line]["critical_value"] = critical_value
        out = open(folder+"run/"+cell_line+"_"+method+".expr", "wb")
        for gene in cell_line_dictionarry[cell_line]:
            if cell_line_dictionarry[cell_line][gene] > cell_line_dictionarry[cell_line]["critical_value"]:
                out.write(gene.replace('"','') +"\n")
        out.close()


cell_line_extractor("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/", "gex.csv", "SD")
