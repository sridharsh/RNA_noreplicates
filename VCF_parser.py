# importing the required modules
from bs4 import BeautifulSoup
import requests
import re
import vcf



# reading in the input VCF file
vcf_reader= vcf_reader = vcf.Reader(open('Challenge_data.vcf', 'r'))
# opening the output file to write the parsed data
f = open('VCF_output.txt', 'a')
# this is the  ExAC browser base url, to web scrape the allele frequencies
base_url = "http://exac.hms.harvard.edu/rest/variant/"

# this is the regex which is used further in the code to parse into the html page and find the allele frequency
pattern = re.compile("allele_freq...(\d\.\d+)")

# this dictionary is used to provide detailed description of the variant types and subtypes
full_forms = {'del': "deletion",
              'ins': "insertion",
              'ts': "transition",
              'tv': "transversion",
              "snp" : "Single Nucleotide Polymorphism ",
              'mnp': "Multi Nucleotide Polymorphism ",
              "unknown": "unknown"}

# writing a header for the output file
header = "Variant type\t " \
         "Variant subtype\t" \
         "Depth of sequence coverage\t " \
         "Number of reads supporting the variant \t " \
         "percent of reads in variant vs reference\t" \
         "ExAC allele-freq\t" \
         "Other information\n"
f.writelines(header)

for record in vcf_reader:
    #print(record)
    variant_type = None  # initializing the variant type in a full form
    for type in record.INFO["TYPE"]:
        if type == "complex":
                variant_type = full_forms[record.var_subtype]
        else:
            variant_type = full_forms[type]
        # print(variant_type)

        for variant in record.ALT:  # parsing through the alternative alleles

            # initializing the ratio between variants and reference reads to null,
            # furthermore, percent is the percentage of the ratio
            # other is the extra information I have provided about the values and variants

            ratio = percent_alt_ref = other = None
            if record.var_type != "indel":

                # First, let's parse though non-indels (SNPs)
                # this part of the code, makes an extension (variant ID),
                # which will be added to the base url to access the ExAC web page,
                # for the web scraping the allele frequency

                ext = (str(record.CHROM) + "-" + str(record.POS) + "-" + str(record.REF) + "-" + str(variant))
                page = requests.get(base_url + ext)
                soup = BeautifulSoup(page.content, 'html.parser')
                string = list(soup.children)[0]

                # using the REGEX pattern to find allele frequency
                match = re.finditer(pattern, string)
                for m in match:
                    allele_freq = (m.group(1))

                    # I have used the try block here, to avoid mathematical error caused
                    # While calculating the percentage of the ratio, there are possibilities of
                    # a zero value for the number of reads supporting the reference allele,
                    # therefore, the except statement feeds in a value "inf" instead of throwing
                    # the zero division error at the user.
                    try:
                        ratio = float(("".join(str(record.INFO['AO']))).strip("[]"))/float("".join(str(record.INFO['RO'])))
                        percent_alt_ref = round(ratio*100, 2)
                        other = ""
                    except ZeroDivisionError: # when the observed count is zero for reference allele
                        percent_alt_ref = "inf"
                        other = "zero reads support ref"

                    f.writelines(variant_type + "\t" + full_forms[record.var_subtype] + "\t" +
                                 str(record.INFO['DP']) + "\t" +
                                 str(("".join(str(record.INFO['AO']))).strip("[]")) + "\t" +
                                 str(percent_alt_ref) +  "\t" +
                                 str(round(float(allele_freq),2)) + "\t" +
                                 str(other) + "\n")

            elif record.var_type == "indel":
                # Now, let's parse though indels (insertions and deletions)

                for ao in record.INFO['AO']:
                    try:
                        ratio = float("".join(str(ao)).strip("[]")) / \
                                float("".join(str(record.INFO['RO'])))
                        percent_alt_ref = round(ratio * 100, 2)
                        other = ""
                    except ZeroDivisionError: # when the observed count is zero for reference allele
                        percent_alt_ref = "inf"
                        other = "zero reads support ref"

                    if record.var_subtype == "ins":
                        # this part of the code, makes an extension (variant ID) only for insertions,
                        # from the bases that are inserted instead of the ref allele
                        # which will be added to the base url to access the ExAC web page,
                        # for the web scraping the allele frequency

                        other = other + "zero reads support ref"

                        for alt in record.ALT:
                            ext = (str(record.CHROM) + "-" + str(record.POS) + "-" + str(record.REF[0]) + "-" +
                                    str(str(alt)[1]))
                            page = requests.get(base_url + ext)
                            soup = BeautifulSoup(page.content, 'html.parser')
                            string = list(soup.children)[0]
                            match = re.finditer(pattern, string)
                            for m in match:
                                allele_freq = (m.group(1))

                                # writing into output file
                                f.writelines(record.var_type + "\t" + variant_type + "\t" +
                                                 str(record.INFO['DP']) + "\t" +
                                                 str(("".join(str(record.INFO['AO']))).strip("[]")) + "\t" +
                                                 str(percent_alt_ref) + "\t" +
                                                 str(round(float(allele_freq), 2)) + "\t" +
                                                 str(other) + "\n")
                    elif record.var_subtype == "del":
                        other = other + "deletion with no ExAC allele freq, so, used AC from VCF file"
                        for alt in record.ALT:
                            f.writelines(str(record.var_type) + "\t" + variant_type + "\t"+ str(record.INFO['DP']) + "\t"+
                                             str(("".join(str(ao))).strip("[]")) + "\t" +
                                             str(percent_alt_ref) + "\t" +
                                             str(round(float(("".join(str(record.INFO["AF"]))).strip("[]")), 2)) + "\t" +
                                             str(other) + "\n")
            else:
                print(record.var_type) # To check if there are any other values

f.close()

