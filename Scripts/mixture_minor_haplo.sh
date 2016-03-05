#!/bin/bash
#
# This script is provided as an example of deriving minor component haplotypes contained in a mixture.
# It is assumed that the normal STR-Seq genotyping/haplotyping pipeline has already been run for
#  both the minor component sample, and the mixture, and intermediate files are available.
#
# Note that running these scripts will create some files which have the same name as those generated
#  in the normal genotyping pipeline since same scripts are run (just with different SNV positions).
#  Therefore you will want to run in a different directory from the one that initial genotyping was run
#  in, and just symlink to the required input files.
#
STRSEQ_DIR='/path/to/strseq'
STR_SCRIPT_DIR="${STRSEQ_DIR}/Scripts/python"
#
fprefix_minor='L002060_S1fixed_L001'
fprefix_mixture='L002061_S4fixed_L001'
#
#Input files needed for this script:
#  <fprefix_minor>_R1.STR_SNV.final.txt
#  <fprefix_minor>_R2.ptag.st.nfilter.trim.st.filter.snv_coords.txt
#  <fprefix_mixture>_R1.STRln_detail.txt
#  <fprefix_mixture>_R2.ptag.st.nfilter.trim.st.bam (+ associated .bai)

mixture_r2_bam=${fprefix_mixture}_R2.ptag.st.nfilter.trim.st.bam
minor_snv_coords=${fprefix_minor}_R2.ptag.st.nfilter.trim.st.filter.snv_coords.txt

# Extract/summarize minor component SNV positions from mixture R2 bam (output: ${fprefix_mixture}_R2.SNV_detail.txt)
python ${STR_SCRIPT_DIR}/pstr_extract_R2snv.py $mixture_r2_bam $minor_snv_coords all

# Merge STR/SNV (output: ${fprefix_mixture}_R1.STR_SNV.detail.txt; ${fprefix_mixture}_R1.STR_SNV.summary.txt)
python ${STR_SCRIPT_DIR}/pstr_merge_str_snv.py ${fprefix_mixture}_R1.STRln_detail.txt ${fprefix_mixture}_R2.SNV_detail.txt all

# Haplotyping of minor component in mixture (output: ${fprefix_mixture}_R1.STR_SNV.minor_haplotypes.txt)
python ${STR_SCRIPT_DIR}/pstr_minor_haplotypes.py ${fprefix_mixture}_R1.STR_SNV.summary.txt ${fprefix_minor}_R1.STR_SNV.final.txt
