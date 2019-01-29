import os, sys
import re
import datetime


def Directory_Maker(sample_file):
    sample_path = "/sc/orga/projects/pacbio/Chromium/"
    samplefile = open(sample_file, "r")
    title = None
    for line in samplefile:
        line = re.sub('[^A-Za-z0-9]+', '\t', line)
        line = line.strip("\n")
        if "CL" in line:
            line = line.split("\t")
            x = re.sub('[^0-9]+', '', line[2])
            p = [line[1], x]
            title = "".join(p)

    now = datetime.datetime.now()
    path = sample_path + str(now.strftime("%Y_%m_%d_" + title))
    if not os.path.exists(path):
        os.makedirs(path)
    os.chdir(path)
    if not os.path.exists("3.0"):
        os.makedirs('3.0')
        os.symlink(os.path.dirname(sample_file), "Raw")
    return ((path+"/3.0/"), path)


def printing_bsubs(job_name, fastq_path, sample_name, fastq_filename, lsf_outfile):

    sys.stdout = lsf_outfile
    print "#BSUB -J", job_name+"_"+sample_name
    print "#BSUB -P acc_apollo"
    print "#BSUB -q premium"
    print "#BSUB -n 4"
    print "#BSUB -W 100:00"
    print "#BSUB -R rusage[mem=12000]"
    print "#BSUB -R span[hosts=1]"
    print "#BSUB -o", job_name+"_"+sample_name+"_log.out"
    print "#BSUB -e", job_name+"_"+sample_name+"_log.err"
    print "\ncd ", fastq_path
    print "module load cellranger/3.0.1"
    print "cellranger mkfastq --run="+fastq_path + "../Raw" \
                       " --csv=sample_sheet.csv --localcores=12 --localmem=47"
    print "cellranger count --id=BC_179_PBMC --sample=" + sample_name + " --fastqs=" + fastq_filename + "/outs/fastq_path " \
        "--transcriptome=/sc/orga/projects/pacbio/Chromium/References/refdata-cellranger-GRCh38-3.0.0 --localcores=12 " \
        "--localmem=47"
    lsf_outfile.close()
    return


def SampleSheet_Maker(sample_file, sample_name, index_file):

    (path, lsf_path) = Directory_Maker(sample_file)
    samplefile = open(sample_file, "r")
    sheet_header = ['Lane', 'Sample', 'Index']
    sheet_values = ['*', sample_name]
    index_file = open(index_file, "r")
    indices = {}
    fetched_idx = []

    outfile = open(path + "sample_sheet.csv", "w")
    lsf_outfile = open(lsf_path + "/run.lsf", "w")

    for indx in index_file:
        indx = indx.strip()
        indx = indx.split(",")
        indices[indx[0]] = indx[1:]
    job_name = fastq_filename = None
    for line in samplefile:
        line = re.sub('[^A-Za-z0-9]+', '\t', line)
        line = line.strip("\n")

        if "CL" in line:
            line = line.split("\t")
            job_name = line[6]
            fastq_filename = line[-8]
            for k, v in indices.items():
                if line[9] in v:
                    fetched_idx.append(k)

    if len(set(fetched_idx)) > 1:
        print("Error! More than one index found!")
    else:
        for x in set(fetched_idx):
            sheet_values.append(x)

    sys.stdout = outfile
    print(",".join(sheet_header))
    print(",".join(sheet_values))
    outfile.close()


    printing_bsubs(job_name, path, sample_name, fastq_filename, lsf_outfile)

    return


if __name__ == "__main__":
    SampleSheet_Maker(sample_file=str(sys.argv[1]), sample_name=str(sys.argv[2]), index_file=str(sys.argv[3]))

