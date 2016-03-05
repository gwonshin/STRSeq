#!/bin/bash
#
RESOURCE_DIR='/path/to/resources'
STRSEQ_DIR='/path/to/strseq'
TMP_DIR='/path/to/tempdir'
JAVA_EXE='/usr/bin/java'

# This path assumes igenomes directory structure, and 1KGenomes fasta file (human_g1k_v37.fasta)
# If your path is different, modify this GENOME_REF_DIR and the --fasta-reference in the freebayes
#   steps below, as needed
GENOME_REF_DIR="${RESOURCE_DIR}/GenomeRef/Homo_sapiens"

BEDTOOLS_DIR="${RESOURCE_DIR}/tools/bedtools2-2.25.0/bin"
VCFTOOLS_DIR="${RESOURCE_DIR}/tools/vcftools_0.1.11/bin"
VCFLIB_DIR="${RESOURCE_DIR}/tools/vcflib/bin"
BAMUTIL_DIR="${RESOURCE_DIR}/tools/bamUtil_1.0.13/bamUtil/bin"
PICARD_DIR="${RESOURCE_DIR}/tools/picard-tools-1.97"
FREEBAYES_EXE="${RESOURCE_DIR}/tools/Freebayes/freebayes"
SAMTOOLS_EXE="${RESOURCE_DIR}/tools/samtools-0.1.18/samtools"

STR_INFO_DIR="${STRSEQ_DIR}/Resources"
STR_SCRIPT_DIR="${STRSEQ_DIR}/Scripts/python"

# Function to print help
print_usage()
{
	   echo "Usage: $(basename "$0") -b {R1bam} -p {OSSeq_pool} [-r {rd_len} -k {flank_size} -v {snv_filter} ]";
	   echo "Where -b R1 bam file must be specified (and must have matching _R2)";
	   echo "      -p pool file must be specified (eg OS0037)";
	   echo "      -r read length, default = 150";
	   echo "      -k flank size, default = 15";
	   echo "      -v R2 snv filter [all|alt], default = alt";     
	   echo "      -d debug (do not delete intermediate files)";
	   return
}

# set default parameters
rdlen=""; flank=""; snvs="";
default_rdlen=150
default_flank=15
default_snvs='alt'
clean='Y'
re_int='^[0-9]+$'

# Parse command line options
OPTIND=1
while getopts "b:p:r:k:v:dh" OPT
do
  case "$OPT" in
    b) r1bam="$OPTARG";;
    p) pool="$OPTARG";;
    r) rdlen="$OPTARG";;
    k) flank="$OPTARG";;
    v) snvs="$OPTARG";;
    d) clean="N";;
    h) print_usage; exit 1;;
   \?) print_usage; exit 1;;
    :) echo "Option -$OPTARG requires an argument."; print_usage; exit 1;;
  esac
done

if [ "$r1bam" == "" ]; then
  echo "R1 bam file must be specified"
  print_usage; 
  exit 1
fi

ext=${r1bam: -3}
if [ "$ext" != 'bam' ]; then
  echo "Input file set must have .bam extension"
  print_usage;
  exit 1
fi

if [ "$pool" == "" ]; then
  echo "OSSeq pool must be specified"
  print_usage; 
  exit 1
fi

if [ "$rdlen" == "" ]; then
  rdlen=$default_rdlen
elif ! [[ $rdlen =~ $re_int ]]; then
  echo "Read length provided ($rdlen) is not an integer"; exit 1
fi

if [ "$flank" == "" ]; then
  flank=$default_flank
elif ! [[ $flank =~ $re_int ]]; then
  echo "Flank size provided ($flank) is not an integer"; exit 1
fi

if [ "$snvs" == "" ]; then
  snvs=$default_snvs
fi

r2bam=${r1bam/_R1/_R2}

r1base=${r1bam%%.*}
r2base=${r2bam%%.*}

export FLANK_SIZE=$flank 

snvs='alt'
haplo='major'
if [ "$snvs" == "all" ]; then
  haplo='all'
fi 

echo "Running STR genotyping with bams: $r1bam, $r2bam"
#R1: Call SNVs in R1 flanking regions
${FREEBAYES_EXE} --pvar 0.05 --no-indels --standard-filters --min-coverage 0 \
    --min-base-quality 20  --min-mapping-quality 30 \
	--fasta-reference ${GENOME_REF_DIR}/1KGenomes/build37/Sequence/BWAIndex/human_g1k_v37.fasta \
    --bam $r1bam --vcf $r1base.vcf

#R1: Extract variants within STR flanking regions
#    Remove any mnp or indel variants called by FreeBayes, and reformat output -> selected columns
${BEDTOOLS_DIR}/intersectBed -a $r1base.vcf \
    -b $STR_INFO_DIR/${pool}.5prflank.st.bed -b $STR_INFO_DIR/${pool}.3prflank.st.bed \
    -names '5pr' '3pr' -wa -wb >$r1base.flank_snv.txt
	
python ${STR_SCRIPT_DIR}/str_flank_alleles.py $r1base.flank_snv.txt

#R1: Genotyping using length-based method
python $STR_SCRIPT_DIR/str_lengths_R1ref.py $STR_INFO_DIR/${pool}.str_probes.txt $STR_INFO_DIR/${pool}.str_info.txt \
    $r1base.flank_alleles.txt $r1bam

python $STR_SCRIPT_DIR/str_ctlen_genotype.py $r1base.STRln_summary.txt $r1base.STRln_probects.txt $STR_INFO_DIR/${pool}.str_info.txt

#R2-SNV Phasing: Using qname, extract R2 mates for R1s which contain STR reads
cut -f5 $r1base.STRln_detail.txt | sort > $r1base.qnames.tmp

${JAVA_EXE} -jar ${PICARD_DIR}/FilterSamReads.jar INPUT=$r2bam TMP_DIR=$TMP_DIR VALIDATION_STRINGENCY=LENIENT \
	    FILTER=includeReadList READ_LIST_FILE=$r1base.qnames.tmp \
	    SORT_ORDER=coordinate CREATE_INDEX=true WRITE_READS_FILES=false \
		OUTPUT=$r2base.nfilter.bam

#R2-SNV Phasing: Mask synthetic DNA positions (40b R2 probes), and last R2 base
#  Last base is prone to false positive variants
${BAMUTIL_DIR}/bam trimBam $r2base.nfilter.bam $r2base.nfilter.trim.bam --left 40 --right 1

#R2-SNV Phasing: Sort/index bam file
${SAMTOOLS_EXE} sort $r2base.nfilter.trim.bam $r2base.nfilter.trim.st
${SAMTOOLS_EXE} index $r2base.nfilter.trim.st.bam

#R2-SNV Phasing: Call SNVs in R2
${FREEBAYES_EXE} --pvar 0.05 --no-mnps --no-complex --min-coverage 3 \
      --min-mapping-quality 25 --min-base-quality 15 --min-supporting-mapping-qsum 90 --min-supporting-allele-qsum 60 \
	  --fasta-reference ${GENOME_REF_DIR}/1KGenomes/build37/Sequence/BWAIndex/human_g1k_v37.fasta \
      --bam $r2base.nfilter.trim.st.bam  --vcf $r2base.nfilter.trim.st.vcf

#R2-SNV Phasing: Extract R2 SNV positions which are in regions defined by bed file (around probe, but not in STR repeat)
#  Additionally apply quality filter using vcffilter
${VCFTOOLS_DIR}/vcftools --vcf $r2base.nfilter.trim.st.vcf \
    --bed ${STR_INFO_DIR}/${pool}_genomic_${rdlen}b.noSTR_plus5b.bed \
    --thin 6 --remove-filtered-all --remove-indels \
    --recode --recode-INFO-all \
    --out $r2base.nfilter.trim.st.filter
	
${VCFLIB_DIR}/vcffilter -f "QUAL > 1 & QUAL / AO > 8" \
    $r2base.nfilter.trim.st.filter.recode.vcf > $r2base.nfilter.trim.st.filter.vcf
rm $r2base.nfilter.trim.st.filter.recode.vcf

#R2-SNV Phasing: Extract pertinent vcf columns and reformat to text
${VCFTOOLS_DIR}/vcftools --vcf $r2base.nfilter.trim.st.filter.vcf --get-INFO TYPE --get-INFO AF \
	--out $r2base.nfilter.trim.st.filter.snv_coords
mv $r2base.nfilter.trim.st.filter.snv_coords.INFO $r2base.nfilter.trim.st.filter.snv_coords.txt

#R2-SNV Phasing: Extract R2 bam reads which cover an SNV position
python $STR_SCRIPT_DIR/pstr_extract_R2snv.py $r2bam $r2base.nfilter.trim.st.filter.snv_coords.txt $snvs
  
#R2-SNV Phasing: Merge R1 STRs with R2 SNVs
python $STR_SCRIPT_DIR/pstr_merge_str_snv.py $r1base.STRln_detail.txt $r2base.SNV_detail.txt 
  
#R2-SNV Phasing: Genotyping
python $STR_SCRIPT_DIR/pstr_genotyping.py $r1base.STR_SNV.summary.txt $STR_INFO_DIR/${pool}.str_info.txt \
    $r1base.STRln_probects.txt 

#Haplotype counts
python $STR_SCRIPT_DIR/pstr_haplotype_cts.py $r1base.STR_SNV.final.txt $haplo

#Cleanup
if [ "clean" == "Y" ]; then
  rm -f $r1base.flank_snv.txt
  rm -f *.nfilter*.bam
  rm -f *.SNV_detail.txt
  rm -f *nfilter.trim*.bam
  rm -f *nfilter*.bai
  rm -f *.tmp
fi
