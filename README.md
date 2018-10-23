# NGS_scripts
Tools for NGS analysis

# VCF parser program

This is a variant annotation tool, to output a table annotating each variant in VCF format file. 

Each variant is annotated with the following pieces of information:
  1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.):
     If there are multiple possibilities, annotate with the most deleterious possibility.
  2. Depth of sequence coverage at the site of variation.
  3. Number of reads supporting the variant.
  4. Percentage of reads supporting the variant versus those supporting reference reads.
  5. Allele frequency of variant from Broad Institute ExAC Project API

     (API documentation is available here: http://exac.hms.harvard.edu/)
  6. Additional optional information from ExAC that you feel might be relevant.

# ChIP overlaps program

This program compares ChIPseq samples of different cell types and identifies overlaps between the MAnorm results or Island lists.

Input: Island list or MAnorm result 
Output: VennDiagram showing the overlapping peaks in the conditions
Note: The program uses euler package to plot the venndiagram, therefore, changes should be made if the user needs to compare >2 conditions or cell types.
