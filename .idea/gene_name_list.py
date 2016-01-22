import string
def search_for_gene_names(startset, infile, sep):
    """
    Makes a set froma  feautre file
    :param: startset: intial set if any
    :param infile: inport feautre file
    :param sep spearator in infile
    :return: set of names
    """
    inp=open(infile)
    if not startset:
        startset=set()
    line = inp.readline()
    line = line.split(sep)
    for unit in line[1:]:
        startset.add(unit)
    return startset


def string_eater(string_file):
    """
    This function eats and digest the string flat file of a particular species and gives the entities in one set.
    :param string_file: the string flat file with multiply weights
    :return: set of identifiers
    """
    inp= open(string_file)
    string_identifiers=set()
    inp.readline()
    for line in inp:
        line=line.split(" ")
        if int(line[7]) >0 or int(line[6])>0:
            a = string.replace(line[0], "9606.", "")
            b = string.replace(line[1], "9606.", "")
            string_identifiers.add(a)
            string_identifiers.add(b)
    return string_identifiers


#nameset=search_for_gene_names(0, "/home/dm729/PycharmProjects/Graph_similarity_measures/cellline_disease_inference_fingerprints.csv", ",")
#nameset=search_for_gene_names(nameset, "/home/dm729/PycharmProjects/Graph_similarity_measures/cell_line_gene_distance_fingerprints.csv", ",")

out= open("/home/dm729/PycharmProjects/Graph_similarity_measures/string_ENSP.txt","wb")

nameset=string_eater("/home/dm729/PycharmProjects/Graph_similarity_measures/9606.protein.links.detailed.v10.txt")
for a in nameset:
    out.write(a)
    out.write("\n")
out.close()


