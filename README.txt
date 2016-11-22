STR-Seq: a massively parallel microsatellite sequencing and genotyping technology

Technology described in publication xyz [include link]
Alignment files are accessible at NCBI Short Read Archive(SRA) with id: SRP071335

Scripts and resource files provided here to rerun any analysis performed for publication

Basic processing can be run in one of two ways.
1/ Run strseq_genotyping.sh with the following parameters:
   -b <R1 bam file>  (Bam file must be indexed, and R2 bam file with same name (except for R2 vs R1) also in same dir and indexed)
   -p <pool>         (Use OS0037 for assay1 and OS0035 for assay2)
   -r <rdlength>     (Alignment read length: 101 if read length is 101, otherwise use 150 which is default)
   -k <flanksize>    (Size of 5' and 3' STR flanking region, max 15.  Use 8 if read length is 101, otherwise 15)
   -v <alt|all>      (Haplotype only Ref/Alt bases at SNV positions (alt), or all bases at SNV positions (all))
  
2/ Run bpipe_str_genotyping.sh (requires installation of bpipe - see below), with the following parameters:
   -f <bam file set, eg *.bam>  (Bam files must be indexed and exist in R1/R2 pairs)
   -p <pool>         (Use OS0037 for assay1 and OS0035 for assay2)
   -r <rdlength>     (Alignment read length: 101 if read length is 101, otherwise use 150 which is default)
   -k <flanksize>    (Size of 5' and 3' STR flanking region, max 15.  Use 8 if read length is 101, otherwise 15)
   -v <alt|all>      (Haplotype only Ref/Alt bases at SNV positions (alt), or all bases at SNV positions (all))
   
Haplotyping of minor component in a mixture is run subsequent to the basic genotyping/haplotyping of both the component (control) 
and the mixture sample, as follows:
3/ Run mixture_minor_haplo.sh  Note this is an example script, it will need to be modified to reflect the actual mixture and minor 
                               component files you wish to run.  More detail provided in comments within the shell script.
							   
Open source tools required (version used in parentheses):
---------------------------------------------------------
BEDTOOLS (v2-2.25.0)
VCFTOOLS (v0.1.11)
VCFLIB
BAMUTIL (v1.0.13)
PICARD  (v1.97)
FREEBAYES (v0.9.21-19)
SAMTOOLS (v0.1.18)

Python dependencies:
--------------------
pysam, distance, numpy, pandas

Shell script variables:
-----------------------
The provided shell scripts will need to be modifed to specify paths to your installation of open source tools, and STR-Seq resources/scripts

Bpipe installation:
-------------------
Overview and download instructions for bpipe are available at: https://github.com/ssadedin/bpipe
Version 0.9.8.7 used in our testing.

