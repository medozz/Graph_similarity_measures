"""
Author: Dezso Modos
This little script will eat the signor whole data tsv file and construct an ncol protein interaction file from it,
which Igraph can consume.
"""
def signor_eater(signor_file, out_file):
    """
    Only the uniprot intearcations are counted . Other moleculaes are not incorporeated.
    :param signor_file: th place of the singor file
    :param out_file: the place and name of the output file
    """
    inp = open(signor_file, "rb")
    edges_set=set() #to delete duplications
    for line in inp:
        line = line.strip()
        line = line .split("\t")
        if line[1]=="PROTEIN" and line[5]=="PROTEIN":
            edge = line[1]+" "+line[5]
            edges_set.add(edge)
    inp.close()

    out=open(out_file,"wb")
    for edge in edges_set:
        out.write(edge+"\n")
    out.close()


signor_eater("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/Signor_2016_02_05.tsv", "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/Signor_2016_02_05.ncol")