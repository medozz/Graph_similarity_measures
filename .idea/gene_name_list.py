import string
def search_for_gene_names(startset, infile, sep):
    """
    Makes a set from a feautre file
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
            #a = string.replace(line[0], "9606.", "")
            #b = string.replace(line[1], "9606.", "")
            a = line[0]
            b = line[1]
            string_identifiers.add(a)
            string_identifiers.add(b)
    return string_identifiers


def outwrite_names(nameset, out_file):
    out=open(out_file,"wb")
    for a in nameset:
        out.write(a)
        out.write("\n")
    out.close()



def string_writer(string_file, uniprot_translation, SwissProtUniProtIds, outfile_up_string_connection, outfile_no_ids):
    uniprot_inp=open(uniprot_translation)
    ENSP_to_UP = {}
    for line in uniprot_inp:
        line=line.split("\t")
        for enspid in line[0].split(","):
            if line[1] in SwissProtUniProtIds:
                if enspid not in ENSP_to_UP:
                    ENSP_to_UP[enspid] = set()
                    ENSP_to_UP[enspid].add(line[1])
                else:
                    ENSP_to_UP[enspid].add(line[1])
    print "ENSP TO UP in memory"
    print len(ENSP_to_UP)
    inp=open(string_file,"rb")
    out=open(outfile_up_string_connection, "wb")
    inp.readline()
    failure=set()
    edges=set()
    for line in inp:
        line=line.split(" ")
        if int(line[7]) >0 or int(line[6])>0:
            a = line[0]
            b = line[1]
            #a = string.replace(line[0], "9606.", "")
            #b = string.replace(line[1], "9606.", "")
            #print ENSP_to_UP[a]
            try:
                for up1 in ENSP_to_UP[a]:
                    try:
                        for up2 in ENSP_to_UP[b]:
                            edges.add(up1+" "+up2)
                    except:
                        failure.add(b)
            except:
                failure.add(a)
    out2=open(outfile_no_ids,"wb") "/home/dm729/PycharmProjects/Graph_similarity_measures/no_UP_STRING.txt","wb")
    for fail in failure:
        out2.write (fail+"\n")
    print failure
    print len(failure)
    for edge in edges:
        out.write(edge+"\n")
    out2.close()
    out.close()

def uniprotin(uniprotfile):
    inp = open(uniprotfile)
    inp.readline()
    SwissProtUniProtIds=set()
    for line in inp:
        line=line.split("\t")
        SwissProtUniProtIds.add(line[0])
    return SwissProtUniProtIds



def reactome_eater(reactomefile, outfile, SwissProtUniProtIds):
    inp = open(reactomefile)
    all_reactions=set()
    inp.readline() #header
    for line in  inp:
        line=line.split("\t")
        a=string.replace(line[0].strip(), "UniProt:","")
        b=string.replace(line[3].strip(), "UniProt:","")
        if a in SwissProtUniProtIds and b in SwissProtUniProtIds:
            interaction=a+" "+b
            all_reactions.add(interaction)
    out= open(outfile,"wb")
    for interaction in all_reactions:
        print interaction
        out.write(interaction)
    out.close()

SwissProtUniProtIds=uniprotin("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/uniprot-homo+sapiens.tab")

reactome_eater("/home/dm729/PycharmProjects/Graph_similarity_measures/homo_sapiens.interactions.txt",
               "/home/dm729/PycharmProjects/Graph_similarity_measures/Reactome_2016_01_22.ncol")

#nameset=search_for_gene_names(0, "/home/dm729/PycharmProjects/Graph_similarity_measures/cellline_disease_inference_fingerprints.csv", ",")
#nameset=search_for_gene_names(nameset, "/home/dm729/PycharmProjects/Graph_similarity_measures/cell_line_gene_distance_fingerprints.csv", ",")

string_writer("/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/homo_sapiens.interactions.txt"
        "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/STRING_ID_TO_UP.tab",
              SwissProtUniProtIds,
              "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/UP_string_only_SP_2016_01_25.ncol",
              "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/no_UP_STRING_SP_2016_01_26.txt")
"""
out_file="/home/dm729/PycharmProjects/Graph_similarity_measures/string_string_ids.txt"

nameset=string_eater("/home/dm729/PycharmProjects/Graph_similarity_measures/9606.protein.links.detailed.v10.txt")

outwrite_names(nameset,out_file)
"""