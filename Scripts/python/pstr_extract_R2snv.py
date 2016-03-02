#!/usr/bin/python

# File: pstr_extract_R2snv.py
# Name: Sue Grimes
# Desc: Script extracts SNV reads from input bam file, using the information
#       contained in the variant .txt file (derived from .vcf)
#
# 6/12/2015: Original version
# 10/6/2015: Modify to use file_open method for input/output files
# 2/8/2016:  Modify to deal with vcf GT values 2/1 or 1/2 etc (multiple alt alleles)
#            Also with case where GT is null (./.)

import os, sys, csv, imp, MySQLdb, pysam, msi_str as msi

script_name = os.path.basename(__file__)
user_home   = os.path.expanduser("~")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if len(sys.argv) < 3:
  print "Usage: ", script_name, "<R2bam_file> <variant_file> [all|alt] [debug]"
  sys.exit(1)

snv_filter = 'alt'
if len(sys.argv) > 3 and sys.argv[3] in ['all', 'both']:
  snv_filter = sys.argv[3]

debug = False
if len(sys.argv) > 4 and sys.argv[4] == 'debug':
  debug = True 
 
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Open/initialize output files and general variables                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
csv.register_dialect('tab_delim', delimiter='\t', doublequote=False, quotechar='', lineterminator='\n', escapechar='', quoting=csv.QUOTE_NONE)

sam_fn = sys.argv[1]
sam_or_bam = sam_fn[-3:]
if os.path.isfile(sam_fn) and os.access(sam_fn, os.R_OK):
  sam_input = pysam.Samfile(sam_fn,'rb') if sam_or_bam == 'bam' else pysam.Samfile(sam_fn,'r')
else:
  print "Unable to open {0} file for input: {1}".format(sam_or_bam, sam_fn)
  sys.exit(1)

snv_input = msi.open_file(sys.argv[2], 'r')

bam_fbase = sam_fn.split('/')[-1].split('.')[0]
r2_output = msi.open_file(bam_fbase + '.SNV_detail.txt', 'w')
snv_out_csv = csv.writer(r2_output, dialect='tab_delim')

print "\n**Running {0}, with R2 bam input: {1}".format(script_name, sam_fn)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# General methods                                                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
def intersect(a1, a2):
  return list(set(a1) & set(a2))
  
def snv_in_region_to_exclude(snv_index, strand, query_length, debug):
  # Exclude any SNVs in synthetic (probe) region, or end base of read
  if strand == 'm' and (snv_index > (query_length - 40) or snv_index == 0):
    exclude_snv = 'Y'
  elif strand == 'p' and (snv_index < 40 or snv_index == query_length-1):
    exclude_snv = 'Y'
  else:
    exclude_snv = 'N'

  if debug: print "SNV index: {0}, R2 strand: {1}, R2 length: {2}".format(snv_index, strand, query_length)
  return exclude_snv
  
def get_allele_freq(snv_row):
  allele_freq = 0
  if len(snv_row) == 6:
    allele_freq = snv_row[5]
  elif len(snv_row) > 8:
    genotype = snv_row[8].split(':')[0]
    allele_freq = gt_to_af(genotype)
  return allele_freq
  
def gt_to_af(genotype):
# Genotypes: 0|1 or 0\1 or 1|1 or 1/2 etc.
  allele1 = genotype[0]  #First char of string
  allele2 = genotype[2]  #Third char of string
  afreq   = '0'
  
  if allele1 == allele2: 
    if allele1 in ['0', '1']:  # 0-homozygous ref; 1-homozygous alt1
      afreq = allele1
    elif allele1 == '2': # 2-homozygous alt2
      afreq = '0,1'
    elif allele1 == '3':
      afreq = '0,0,1'    # 3-homozygous alt3
    
  elif allele1 in ['0','1'] and allele2 in ['0','1']: # heterozygous ref/alt
    afreq = '0.5'
  else:
    afreq = '0.5,0.5'    # 1/2 or 2/1
  
  return afreq

def ref_or_alt(ref_alt_bases, snv_filter):
  if snv_filter == 'both':
    return [ref_alt_base[0] for ref_alt_base in ref_alt_bases]
  elif snv_filter == 'alt':
    return [ref_alt_base[0] for ref_alt_base in ref_alt_bases if ref_alt_base[1] > 0]
  else: #snv_filter == 'all'
    return ['A', 'C', 'G', 'T']
 
def print_variant(qname, probe_nr, snv_chr, snv_pos, snv_base, ref_bases, probe_region_flag):
  refbases    = [ref_base[0] for ref_base in ref_bases]
  is_or_isnot = 'is not' if probe_region_flag == 'N' else 'IS'
  msg_prefix  = ''       if probe_region_flag == 'N' else '**'
  print "{7} {0}Base: {1} at pos: {2}:{3} {4} in probe {5} synthetic region.  Ref/Alt: {6}".format(msg_prefix, snv_base, snv_chr, snv_pos, is_or_isnot, probe_nr, refbases, qname)
          

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Main program logic                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#-----------------------------------------------------------------------------#
# Read SNVs from input file                                                   #
#-----------------------------------------------------------------------------#
# snv_input (FreeBayes) example:
# CHROM	POS	REF	ALT	TYPE	AF
#1	7142604	C	G	snp	1
#1	10987656	A	C	snp	1
#1	64304937	C	T	snp	0.5
#1	74184409	C	T	snp	1
#1	74229419	G	A	snp	0.5

# snv_input (1000kG) example:
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO							FORMAT		NA12878
#1	52066	rs28402963	T	C	.	PASS	AA=C;DP=84;GP=1:62203;BN=125	GT:GQ:DP	1/0:70:36
#1	695745	rs72631875	G	A	.	PASS	AA=.;DP=124;GP=1:705882;BN=130	GT:GQ:DP	1|0:99:56
#1	742429	rs3094315	G	A	.	PASS	AA=g;DP=132;HM2;GP=1:752566;BN=103	GT:GQ:DP	1|1:99:44
#1	742584	rs3131972	A	G	.	PASS	AA=a;DP=160;HM3;GP=1:752721;BN=103	GT:GQ:DP	1|1:99:60

snv_info = {}
snv_reader = csv.reader(snv_input, dialect='tab_delim')
next(snv_reader)  #Skip header

snv_array = list(snv_reader)
for snv_row in snv_array:
  snv_chr = snv_row[0]
  snv_pos = int(snv_row[1]) # Convert to 0-based positions for matching with pysam/bam  
  ref_bases = snv_row[2]
  alt_bases = snv_row[3]
  allele_freq = get_allele_freq(snv_row)
  if allele_freq not in ['0', '0.5', '1']:
    continue
  ref_alt_bases = msi.snv_ref_alt(ref_bases, alt_bases, allele_freq)
  # Output from msi.snv_ref_alt is array of tuples:  
  # Eg: [('G',0),('T',0.5),('C',0.5)], where first tuple is reference and can be 0 freq.
  if snv_chr in snv_info.keys():
    snv_info[snv_chr][snv_pos] = ref_alt_bases
  else:
    snv_info[snv_chr] = {snv_pos: ref_alt_bases}

#-----------------------------------------------------------------------------#
# Read bam file, and check for SNV from input list                            #
#-----------------------------------------------------------------------------#
sam_rows = sam_input.fetch()
snv_out_csv.writerow(['qname', 'probe_nr', 'ref_chr', 'snv_pos', 'snv_base', 'ref/alt'])

for sam_row in sam_rows:
  sam_positions = sam_row.positions
  sam_chr = '*' if sam_row.is_unmapped else sam_input.getrname(sam_row.reference_id)
  sam_zp_tag = [ tag[1] for tag in sam_row.tags if tag[0] == 'ZP' ]
  probe_nr = sam_zp_tag[0]
  r2_strand = 'm' if sam_row.is_reverse else 'p'
  r2_length = sam_row.query_length
  
  if sam_chr in snv_info.keys():
    #if debug: print "Chr: {0} Reference SNVs: {1}".format(sam_chr, snv_info[sam_chr])
    aligned_snv_pos = intersect(snv_info[sam_chr].keys(), sam_positions)

    if len(aligned_snv_pos) > 0:
      # Get_aligned_pairs() returns position in aligned read, and chromosome coordinate for each aligned base
      # Eg: [(0, 7142487), (1, 7142488), (2, 7142489), (3, 7142490), ...]
      # **NOTE**: These positions are 0-based, but positions from vcf file are 1-based
      aligned_coords = sam_row.get_aligned_pairs()
      if debug: print "\nFound variant pos in R2: {0}, {1}:{2}".format(sam_row.qname, sam_chr, aligned_snv_pos)

      for snv_pos in aligned_snv_pos:
        snv_indices = [seq_idx for (seq_idx, chr_coord) in aligned_coords if chr_coord != None and snv_pos == chr_coord+1]
        if debug: print "SNV indices: {0}".format(snv_indices)
        if debug: print sam_row.qname, sam_chr, sam_positions
        if debug: print "Aligned Coords:{0}".format(aligned_coords)
       
        if len(snv_indices) == 0 or snv_indices[0] == None:
          continue
 
        snv_index = snv_indices[0] #There will only be one value in array due to snv_pos == chr_coord+1 if test
        probe_region_flag = snv_in_region_to_exclude(snv_index, r2_strand, r2_length, debug)
        sam_base = sam_row.query_alignment_sequence[snv_index]
        ref_alt_bases = snv_info[sam_chr][snv_pos]

        if debug: print_variant(sam_row.qname, probe_nr, sam_chr, snv_pos, sam_base, ref_alt_bases, probe_region_flag)
        if probe_region_flag == 'N' and sam_base in ref_or_alt(ref_alt_bases, snv_filter):
          snv_out_csv.writerow([sam_row.qname, probe_nr, sam_chr, snv_pos, sam_base, '/'.join(ref_alt_base[0] for ref_alt_base in ref_alt_bases)])
        
sam_input.close()
snv_input.close()
r2_output.close()