import sys


#run combine_pathway.R before running this python script
def foldchange(gsea_preranked_list = "DEgenes.csv"):

    FC_dict = {}

    file_FC = open(gsea_preranked_list, "r")
    for line1 in file_FC:
        FC_dict[line1.strip().split(",")[1]] = (line1.strip().split(",")[0])

    return FC_dict


def identify_regulation(pathway_counts="pathway_count.txt",
                        outfile="output_updownpathways.tsv"):
    foldchanges_dict = foldchange()
    file_pathway = open(pathway_counts, "r")
    final_output = {}

    for line in file_pathway:

        line = line.strip().split("\t")
        geneid = line[0]
        pathway = line[3]
        if pathway in final_output:
            if geneid in foldchanges_dict:
                if float(foldchanges_dict[geneid]) < 0:
                    final_output[pathway]["downregulated"] += 1
                elif float(foldchanges_dict[geneid]) > 0:
                    final_output[pathway]["upregulated"] += 1
                else:
                    pass
            else:
                pass
        else:
            final_output[pathway] = {"downregulated": 0,
                                     "upregulated": 0}

    outfile = open(outfile, "w")
    sys.stdout = outfile

    print("Pathway","\t", "downregulated", "\t", "upregulated")
    for keys,dict_vals in final_output.items():
        print(keys, "\t", dict_vals["downregulated"], "\t", dict_vals["upregulated"])

    return outfile.close()


identify_regulation()

