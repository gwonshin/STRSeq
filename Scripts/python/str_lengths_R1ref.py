#!/usr/bin/python

# File: str_lengths_R1ref.py
# Name: Sue Grimes
# Desc: Script extracts STR reads from input sam file, using the information
#       contained in the probe and str files (also input) for the relationship 
#       between probe# (tagged in sam file) and str name/motif etc.
#
# 6/15/2015: Original version 
# 6/17/2015: Modify to only select STRs which have the reference flanking seq
# 10/6/2015: Modify to use msi.open_file method for input/output files
# 10/8/2015: Modify to be consistent with str_counts_R1ref.py and to output summary file,
#              but not output the final (matrix) file
# 10/14/2015: Modify to read .bam or .sam file

import os, sys, re, csv, imp, pysam, distance, numpy as np, msi_str as msi

script_name = os.path.basename(__file__)
user_home   = os.path.expanduser("~")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if len(sys.argv) < 4:
  print "Usage: ", script_name, "<probe_info> <str_info> <flank_snvs> <bam_file>"
  sys.exit(1)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Open/initialize output files and general variables                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
csv.register_dialect('tab_delim', delimiter='\t', doublequote=False, quotechar='', lineterminator='\n', escapechar='', quoting=csv.QUOTE_NONE)

probe_input  = msi.open_file(sys.argv[1], 'r')
probe_csv = csv.DictReader(probe_input, dialect='tab_delim')

str_input    = msi.open_file(sys.argv[2], 'r')
str_csv = csv.DictReader(str_input, dialect='tab_delim')

fsnv_input   = msi.open_file(sys.argv[3], 'r')
fsnv_csv = csv.DictReader(fsnv_input, dialect='tab_delim')

sam_fn = sys.argv[4]
sam_or_bam = sam_fn[-3:]
if os.path.isfile(sam_fn) and os.access(sam_fn, os.R_OK):
  sam_input = pysam.Samfile(sam_fn,'rb') if sam_or_bam == 'bam' else pysam.Samfile(sam_fn,'r')
else:
  print "Unable to open {0} file for input: {1}".format(sam_or_bam, sam_fn)
  sys.exit(1)

sam_base     = sam_fn.split('/')[-1].split('.')[0]
r1dtl_output = msi.open_file(sam_base + '.STRln_detail.txt', 'w')
summ_output = msi.open_file(sam_base + '.STRln_summary.txt', 'w')
probe_output = msi.open_file(sam_base + '.STRln_probects.txt', 'w')

r1dtl_hdgs = ['str_name', 'motif', 'probe_nr', 'strand', 'qname', 'chrom', 'pos', 'cigar', 'rlen', 'motif_start', 'motif_cnt',
                           'str_len', 'str_sequence', 'bases_5pr', 'bases_3pr', 'trunc_flag', 'seq', 'mapq', 'qqual']
r1dtl_csv = csv.DictWriter(r1dtl_output, r1dtl_hdgs, dialect='tab_delim')
r1dtl_csv.writeheader()

FLANK_SIZE = msi.FLANK_SIZE
ALLELE2_MIN_PCT = msi.ALLELE2_MIN_PCT

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Main program logic                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#-----------------------------------------------------------------------------#
# Create probe_info and str_info dictionaries, from input csv files           #
# Add alternate flanking sequences if variants in flank positions             #
#-----------------------------------------------------------------------------#
str_info   = msi.dict_from_csv(str_csv, 'str_info')
#print sorted(str_info)

probe_info = msi.dict_from_csv(probe_csv, 'probe_info')

# Modify str_info dictionary to put flanking sequences into array to allow appending later
for str_name, str_vals in str_info.items():
  str_vals['flanking_5pr'] = [str_vals['flanking_5pr']]
  str_vals['flanking_3pr'] = [str_vals['flanking_3pr']]

for frow in fsnv_csv:
  str_name = frow['STRName']
  flank5or3 = 'flanking_' + frow['5or3pr']
  
  snv_index   = int(frow['SNVPos']) - (int(frow['FlankStart'])+1)
  ref_base    = frow['Ref']
  alt_alleles = frow['Alt']
  allele_freq = frow['AF']
  
  bases_wfreq = msi.snv_ref_alt(ref_base, alt_alleles, allele_freq)
  print "Bases w/freq: {0}".format(bases_wfreq)
  
  # Add alternate flanking, for non-reference alleles with frequency > 0
  for (snv_base, afreq) in bases_wfreq[1:]:
    if afreq == 0:
      continue 
    
    # Modify reference and current alt flanking sequences with alternate SNV base
    # Result is that reference, individual mutations and combination mutations are all tested for in flanking sequence
    nr_flanking = len(str_info[str_name][flank5or3])
    for i in range(nr_flanking):
      flank_seq = str_info[str_name][flank5or3][i]
      alt_flank = flank_seq[0:snv_index] + snv_base + flank_seq[snv_index+1:]
      if alt_flank not in str_info[str_name][flank5or3]:
        str_info[str_name][flank5or3].append(alt_flank)
        print "Alt {0} flanking for STR: {1} is: {2}".format(frow['5or3pr'], str_name, alt_flank)
    
#for str_name, str_vals in str_info.items():
#  if len(str_vals['flanking_5pr']) > 1 or len(str_vals['flanking_3pr']) > 1:
#    print str_name, str_info[str_name]
   
#-----------------------------------------------------------------------------#
# Read sam file, and check for motif if probe# is in the input list           #
#-----------------------------------------------------------------------------#
tot_reads = {}; ref_flank_rds = {};
sam_rows = sam_input.fetch()

for sam_row in sam_rows:
  sam_zp_tag = [ tag[1] for tag in sam_row.tags if tag[0] == 'ZP' ]
  probe_nr = sam_zp_tag[0]
  
  if probe_nr in probe_info.keys():
    probe_nr_info = probe_info[probe_nr]
    str_name      = probe_nr_info['str_name']
    str_nm_info   = str_info[str_name]
    strand        = probe_nr_info['strand']
    multi_motif  = ''.join([str_nm_info['str_motif']*str_nm_info['min_repeats']])
    
    probe_info[probe_nr]['probe_reads'] += 1
    msi.increment_dict_ct(tot_reads, (str_name, strand))

    seq_5prto3pr = msi.rev_complement(sam_row.seq) if (probe_nr_info['strand'] == 'M' and sam_row.is_unmapped == True) else sam_row.seq
    found_motif  = msi.motif_search(seq_5prto3pr, sam_row.rlen, str_nm_info['str_motif'], str_nm_info['min_repeats'])

    if found_motif['motif_cnt'] > 0:
      probe_info[probe_nr]['motif_reads'] += 1 
    
    flanking_match = None; flanking_5pr = None; flanking_3pr = None;
    for flanking_5pr in str_nm_info['flanking_5pr']:
      if flanking_match:
        break  #Need to break out of both for loops when match is found

      for flanking_3pr in str_nm_info['flanking_3pr']:
        flanking_match  = re.search(flanking_5pr + r'(.*' + multi_motif + r'.*)' + flanking_3pr, seq_5prto3pr)
        if flanking_match:
          break

    if flanking_match == None:
      continue

    #-----------------------------------------------------------------------------#
    # Write detailed output: All reads with full STR and ref/alt flanking seqs    #
    #-----------------------------------------------------------------------------#
    probe_info[probe_nr]['full_str_rds'] += 1
    #print sam_row.qname, flanking_match.group(1), found_motif
    
    chr_ref = '*' if sam_row.is_unmapped else sam_input.getrname(sam_row.tid)    
    motif_start = msi.chr_pos_motif(sam_row.pos, found_motif['motif_idx'], probe_nr_info['strand'])
    str_sequence = flanking_match.group(1)
    str_len      = len(str_sequence)
    motif_cnt    = "{0:.5g}".format(str_len * 1.0 / len(str_nm_info['str_motif']))

    str_dtl_row = {'str_name': str_name, 'motif': str_nm_info['str_motif'], 'probe_nr': probe_nr, 'strand': probe_nr_info['strand'],
                   'qname': sam_row.qname, 'chrom': chr_ref, 'pos': sam_row.pos, 'cigar': 'NA', 'rlen': sam_row.rlen, 
                   'motif_start': motif_start, 'str_len': str_len, 'motif_cnt': motif_cnt, 'str_sequence': str_sequence, 
                   'bases_5pr': flanking_5pr, 'bases_3pr': flanking_3pr, 'trunc_flag': found_motif['trunc_flag'], 
                   'seq': seq_5prto3pr, 'mapq': sam_row.mapq, 'qqual': sam_row.qqual }  
    r1dtl_csv.writerow(str_dtl_row)

    #-----------------------------------------------------------------------------#
    # Accumulate summary counts (use ref flanking, not actual flanking in summary)#
    #-----------------------------------------------------------------------------#
    key_5pr_3pr = (str_name, strand, str_nm_info['str_motif'], str_nm_info['flanking_5pr'][0], str_nm_info['flanking_3pr'][0])
    if not key_5pr_3pr in ref_flank_rds.keys():
      ref_flank_rds[key_5pr_3pr] = {}    

    msi.increment_dict_ct(ref_flank_rds[key_5pr_3pr], len(str_sequence))
    
#-----------------------------------------------------------------------------#
# Write summary output: All reads with at least minimum STR repeats, between  #
#   exact match to flanking sequences (either reference seq or alt allele)    #
#-----------------------------------------------------------------------------#
str_out_csv = csv.writer(summ_output, dialect='tab_delim')
str_out_csv.writerow(['STR Name', 'Strand', 'Motif', 'Min Rpts', 'Probe Rds', '5pr Flank', '3pr Flank', 'STR Len', 'Motif#', 'STR Rds'])
  
if len(ref_flank_rds) > 0:
  for str_flank_key in sorted(ref_flank_rds):
      str_len_cts = ref_flank_rds[str_flank_key]

      if len(str_len_cts) > 0:
        for str_len_key in sorted(str_len_cts):

          str_name  = str_flank_key[0]
          strand    = str_flank_key[1]
          str_motif = str_flank_key[2]
          flank_5pr = str_flank_key[3]
          flank_3pr = str_flank_key[4]
          str_len   = str_len_key
          motif_rd_cnt  = str_len_cts[str_len_key]
          nr_motifs = "{0:.5g}".format(str_len * 1.0 / len(str_motif))
          str_out_csv.writerow([str_name, strand, str_motif, str_info[str_name]['min_repeats'], tot_reads[(str_name, strand)], 
                                flank_5pr, flank_3pr, str_len, nr_motifs, motif_rd_cnt])
else:
  str_out_csv.writerow(['NA', 'NA', 'NA', 0, 0, 'NA', 'NA', 0, 0, 0])

#-----------------------------------------------------------------------------#
# Write probe counts: Total probe reads per STR and reads which include at    #
#  least the minimum number of repeats of STR motif (irrespective of flanking)#
#-----------------------------------------------------------------------------#
# {'probe_start_pos': 12009236, 'chromosome': '3', 'str_name': 'trf533283', 'strand': 'P', 'probe_reads': 2355, 'probe_nr': 499,  
#   'motif_reads': 419, 'full_str_rds': 156}
probe_hdgs = ['probe_nr', 'chromosome', 'probe_start_pos', 'str_name', 'strand', 'probe_reads', 'motif_reads', 'full_str_rds']
probeout_csv = csv.DictWriter(probe_output, probe_hdgs, dialect='tab_delim')
probeout_csv.writeheader()

for probe_nr in probe_info:
  probe_info[probe_nr]['probe_nr'] = probe_nr
  probeout_csv.writerow(probe_info[probe_nr])
probe_output.close()

r1dtl_output.close()
summ_output.close()

