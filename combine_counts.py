import sys, os


def rpkm_counts(dict_cond, total_counts):

    # RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )

    rpkm_dict = {}

    for k in dict_cond.keys():
        length = dict_cond[k][0]
        counts = dict_cond[k][1]
        rpkm_dict[k] = (counts, counts/((length/1000)*(total_counts/1000000)))

    return rpkm_dict


def merge_file_counts(list_of_files):

    # merging the counts columns based on the geneIDs

    numb_files = len(list_of_files)

    wt_dict = {}
    ko_dict = {}
    total_counts1 = 0
    total_counts2 = 0

    for i in range(0, numb_files - 1):
        file_1 = list_of_files[i]
        file_2 = list_of_files[i + 1]
        header1 = file_1.readline()
        h1 = file_1.readline()
        header2 = file_2.readline()
        h2 = file_2.readline()

        for line1, line2 in zip(file_1, file_2):
            #print(line1)
            #print(line1.strip().split()[-1])
            line1 = line1.strip().split()
            line2 = line2.strip().split()
            #print(line1)
            wt_dict[line1[0]] = (int(line1[-2]), int(line1[-1]))
            ko_dict[line2[0]] = (int(line2[-2]), int(line2[-1]))
            total_counts1 += int(line1[-1])
            total_counts2 += int(line2[-1])
    WT = rpkm_counts(wt_dict, total_counts1)
    KO = rpkm_counts(ko_dict, total_counts2)

    return (WT, KO)

if __name__ == "__main__":

    ADAMTSL2 = open("KO.exon.geneID.txt","r")
    Control = open("WT.exon.geneID.txt","r")
    (WT, KO) = merge_file_counts([Control, ADAMTSL2])

    # printing the results into an output file

    merged_list = []

    outfile = open("rpkm.csv", "w")
    sys.stdout = outfile
    print("GeneID, WT_RPKM, KO_RPKM, WT_counts, KO_counts")
    for wt_k in WT.keys():
        if wt_k in KO.keys():
            merged_list.append(wt_k)
            print(wt_k.upper() + "," +
                  str(WT[wt_k][1]) + "," +
                  str(KO[wt_k][1])+ "," +
                  str(WT[wt_k][0]) + "," +
                  str(KO[wt_k][0]))
        else:
            pass

    #print(merged_list)
    outfile.close()
