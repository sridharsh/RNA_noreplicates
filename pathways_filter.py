from pprint import pprint

# this program fetches the common pathways from two biological conditions given by the user
# the input files are the GO enrichment results from ToppFunn
# the file will be tab-delimited
# make sure you make a list of the file names in the Main section


class common_pathways(object):

    def __init__(self, file1, file2, file_out):

        self.file_1 = open(file1, "r")
        header_1 = self.file_1.readline()

        self.file_2 = open(file2, "r")
        header_2 = self.file_1.readline()

        self.outfile = open(file_out, "w")

        self.list_pathways = {}
        self.count = 0

    def filtering(self, cutoff=None):
        for line in self.file_1:
            line = line.strip()
            line = line.split("\t")
            # print(line[0])
            if cutoff:
                if float(line[7]) <= cutoff:
                    self.list_pathways[line[0]] = 1
                    # print(line[7], line[0])

                else:
                    pass
            else:
                self.list_pathways[line[0]] = 1

        for line in self.file_2:
            line = line.strip()
            line = line.split("\t")
            # print(line[0])
            if line[0] in self.list_pathways.keys():
                self.list_pathways[line[0]] += 1
            else:
                self.list_pathways[line[0]] = 1

        # pprint(self.list_pathways)

        for k,v in sorted(self.list_pathways.items()):
            if v >= 2:
                self.count += 1
                self.outfile.writelines(k+"\n")

        print("Number of common pathways in the given data files: ", self.count)
        print("Total number of pathways in both the files combined: ", len(self.list_pathways))


def main():

    files_HT = ["HT_c2_neg.xls", 
                "HT_c2_pos.xls",
                "HT_c3_neg.xls",
                "HT_c3_pos.xls",
                "HT_h_neg.xls",
                "HT_h_pos.xls", 
                "HT_c5_neg.xls", 
                "HT_c5_pos.xls"]
    files_NT = ["NT_c2_neg.xls",
                "NT_c2_pos.xls",
                "NT_c3_neg.xls",
                "NT_c3_pos.xls",
                "NT_h_neg.xls",
                "NT_h_pos.xls", 
                "NT_c5_neg.xls", 
                "NT_c5_pos.xls"]

    for i in range(len(files_HT)):

        print(files_HT[i], files_NT[i])
        data_load = common_pathways(files_HT[i], 
                                    files_NT[i], 
                                    files_HT[i].replace(".xls", "_")+files_NT[i].replace(".xls", ".txt"))
        data_load.filtering(cutoff=0.01)
        # print("Generating the list of common pathways... \n")


if __name__ == '__main__':
    main()
