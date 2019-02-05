import os, sys
import re
import datetime


def index_dict(index_file):

    '''
    Method to fetch out indices for the given pattern

    Input: index file

    Output: index dictionary with; keys = pattern
                                   values = index
    '''

    indices = {}
    for index in open(index_file, "r"):
        index = index.strip()
        index = index.split(",")
        for i in index[1:]:
            indices[i] = index[0]
    return indices


def printing_bsubs(job_name, fastq_path, sample_name, fastq_filename, lsf_path, i, mode):

    '''
    Method to print out lsf files for the sample file given

       Input: job name
              path to output
              collaborator ID
              name of the folder where the fastq files would be saved
              path to output

       Output:
             incase of multiple samples combined, it will create multiple lsf files; else one single run.lsf file
    '''

    if mode == "single":
        outfile = open(lsf_path + "/run.lsf", "w")
        sys.stdout = outfile
    if mode == "multiple":
        outfile = open(lsf_path + "/run.lsf" + i, "w")
        sys.stdout = outfile

    print ("#BSUB -J", job_name + "_" + sample_name +
           "\n#BSUB -P acc_apollo\n#BSUB -q premium\n#BSUB -n 4"
           "\n#BSUB -W 100:00\n#BSUB -R rusage[mem=12000]"
           "\n#BSUB -R span[hosts=1]\n#BSUB -o" + job_name + "_" + sample_name + "_log.out"
           "\n#BSUB -e", job_name + "_" + sample_name + "_log.err\n")

    print("cd ", fastq_path)

    print("module load cellranger/3.0.1")

    print("cellranger mkfastq --run=" + fastq_path + "../Raw" +
          " --csv=sample_sheet.csv --localcores=12 --localmem=47")

    print("cellranger count "
          "--id=" + sample_name +
          " --sample=" + sample_name +
          " --fastqs=" + fastq_filename +
          "/outs/fastq_path" 
          " --transcriptome=/sc/orga/projects/pacbio/Chromium/References/refdata-cellranger-GRCh38-3.0.0"
          " --localcores=12 --localmem=47")

    outfile.close()

    return


def Directory_Maker(name, outdir="/Users/sridhs01/Documents/"):


    '''
    Method to make directories with date_CL_id as name in Chromium folder and create Raw softlink to the sample data
    and /3.0 folder for the sample sheet

    Input: List of CL ids
           Output directory; which is defaulted to /pacbio/Chromium

    Output: Directories and softlink to the Raw data
    '''

    now = datetime.datetime.now()
    path = outdir + str(now.strftime("%Y_%m_%d_" + "CL_" + "_".join(name)))
    if not os.path.exists(path):
        os.makedirs(path)
    os.chdir(path)
    if not os.path.exists("3.0"):
        os.makedirs('3.0')
        os.symlink(os.path.dirname(outdir), "Raw")
    return path+"/3.0/", path


def sample_sheet(summary, sample_results):

    '''
    Method to produce sample sheet and run.lsf files

    Input: summary list of parsed regex patterns
           sample_results of parsed CL id and indices

    Output: sample_sheet.csv and run.lsf (multiple lsf files if multiple samples are combined)
    '''

    cl_list = []
    for x in sample_results:
        cl_list.append(x[1].strip("CL"))
    path, lsf_path = Directory_Maker(cl_list)
    os.chdir(path)
    outfile = open("sample_sheet.csv", "w")
    sys.stdout = outfile
    print(",".join(('Lane', 'Sample', 'Index')))
    print(summary)
    for results in sample_results:
        print(",".join(results))
    outfile.close()

    if len(summary) > 1:
        for i, value in enumerate(summary):
            printing_bsubs(value[1], path, value[0], value[2], lsf_path,  str(i), mode="multiple")
    else:
        for value in summary:
            printing_bsubs(value[1], path, value[0], value[2], lsf_path, str(0), mode="single")

    return


def regex_parser(sample_file, index_file):

    '''
    Method to parse for regex patterns

    Input: sample file
           index file

    Output: produces sample sheet using the sample_sheet method above
    '''

    re_CL_id = re.compile("(CL).(\d+)")
    re_nt = re.compile("[ATGC]{8}")
    re_collab_id = re.compile("^.*(TD.{5})")
    re_id = re.compile("Homo_Sapien,(.{9})")

    indices = index_dict(index_file)
    summary = []
    sample_results = []

    for line in open(sample_file, "r"):
        Title_found = re.search(re_id, line)
        Collab_found = re.search(re_collab_id, line)
        Index_found = re.search(re_nt, line)
        CL_found = re.search(re_CL_id, line)

        try:
            cl_id = (CL_found.group(1) + CL_found.group(2))
            index = indices[str(Index_found.group(0))]
            collab_id = (Collab_found.group(1))
            title = (Title_found.group(1))
            sample_results.append(("*", cl_id, index))
            summary.append((cl_id, collab_id, title))
        except:
            pass

    return sample_sheet(set(summary), set(sample_results))


if __name__ == "__main__":
    regex_parser(sample_file=str(sys.argv[1]), index_file=str(sys.argv[2]))
    
    
