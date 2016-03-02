#!/usr/bin/python

# File: pstr_minor_haplotypes.py
# Name: Sue Grimes
# Desc: Script takes summary STR/SNV information and a list of haplotypes
#         and outputs haplotypes from the list which occur in the STR/SNV file
#
# 1/20/2016: Original version

import os, sys, csv, imp, pysam, numpy as np, pandas as pd, msi_str as msi
from decimal import Decimal
from collections import OrderedDict

script_name = os.path.basename(__file__)
user_home   = os.path.expanduser("~")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if len(sys.argv) < 3:
  print "Usage: ", script_name, "<str_snv_summary(mix)> <str_snv_final(minor)> [debug]"
  sys.exit(1)

debug = False
if len(sys.argv) > 3 and sys.argv[3] == 'debug':
  debug = True 

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Open/initialize output files and general variables                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
csv.register_dialect('tab_delim', delimiter='\t', doublequote=False, quotechar='', lineterminator='\n', quoting=csv.QUOTE_NONE)

summ_input   = msi.open_file(sys.argv[1], 'r')
minor_input  = msi.open_file(sys.argv[2], 'r')

summ_fbase = sys.argv[1].split('/')[-1].split('.')[0] 

haplo_output = msi.open_file(summ_fbase + '.STR_SNV.minor_haplotypes.txt', 'w')
haplo_csv = csv.writer(haplo_output, dialect='tab_delim')

print "\n**Running {0}, with STR/SNV inputs: {1}, {2}".format(script_name, sys.argv[1], sys.argv[2])

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# General methods                                                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
def extract_haplotypes(snv_base, str_counts):
  # snv_base: 'C'  str_counts: [[10,11], [38,6], 63] 
  # returned tuples: [('C',10,38,0.6,48),('C',11,6,0.1,48)]
  [motif_rpts, snv_reads, tot_rds] = str_counts
  haplo_tuples = []
  
  for i, motif_rpt in enumerate(motif_rpts):
    snv_rd_pct = "{0:.2f}".format((snv_reads[i]*100.0/tot_rds))
    haplo_tuples.append((snv_base, motif_rpt, snv_reads[i], snv_rd_pct, tot_rds))

  if debug: print "Haplotype tuples: {0}".format(haplo_tuples)
  return haplo_tuples
  
def format_haplo_csv(str_snv_key, snv_allele_rds):
  (str_name, snv_chr, snv_pos) = str_snv_key
  alleles_rds_srted = sorted(snv_allele_rds, key=lambda x: int(x[2]), reverse=True) 
  print str_snv_key, alleles_rds_srted
  
  snv_reads = alleles_rds_srted[0][4]  #total reads is repeated for each tuple, so just pick first one

  snv_bases   = [allele_rd[0] for allele_rd in alleles_rds_srted]
  snv_bases_uniq = list(OrderedDict.fromkeys(snv_bases))
  
  str_alleles = [allele_rd[1] for allele_rd in alleles_rds_srted]
  str_alleles_uniq = list(OrderedDict.fromkeys(str_alleles))
  
  haplotypes  = [('-').join([allele_rd[0], allele_rd[1]]) for allele_rd in alleles_rds_srted]

  haplo_rds   = [allele_rd[2] for allele_rd in alleles_rds_srted]
  hap_rd_pcts = [allele_rd[3] for allele_rd in alleles_rds_srted]
  
  return [str_name, snv_chr, snv_pos, snv_reads, ','.join(snv_bases_uniq), ','.join(str_alleles_uniq),
            ','.join(haplotypes), ','.join([str(rds) for rds in haplo_rds]),
            ','.join([hap_rd_pct for hap_rd_pct in hap_rd_pcts])]


#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Main program logic                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#-----------------------------------------------------------------------------#
# Read minor final summary file, and load into dictionary                     #
#-----------------------------------------------------------------------------#
# Input example
#STR name	Motif	Ref Chr	SNV Pos	SNV Allele	SNV Reads	Motif Rpts			Motif Rpt Rds	STR Allele(s)
#D16S539	GATA	16		86386213	A		130			7, 8, 9, 10, 11		1, 1, 1, 8, 119		11
#D16S539	GATA	16		86386213	C		101			12, 13				100, 1				12
#D20S1082	ATA		20		53865850	G		25			11					25					11
#trf281749	ATT		15		92517426	T		58			12, 15, 16, 17, 18	1, 4, 27, 24, 2		16, 17

minor_csv = csv.DictReader(minor_input, dialect='tab_delim')
minor_str = {}

for mrow in minor_csv:
  mhap_key = (mrow['STR name'], mrow['Ref Chr'], int(mrow['SNV Pos']), mrow['SNV Allele'])
  str_alleles = [str_allele for str_allele in mrow['STR Allele(s)'].split(', ')]
  if mhap_key in minor_str:
    minor_str[mhap_key] += str_alleles
  else:
    minor_str[mhap_key] = str_alleles
  
if True: print "\nMinor component STR dictionary:\n{0}".format(minor_str)

#-----------------------------------------------------------------------------#
# Read STR_SNV summary csv file, and load into pandas DataFrame               #
#-----------------------------------------------------------------------------#
# Input example
#str_name	strand	bases_5pr	bases_3pr	motif_cnt	str_ct	ref_chr	snv_pos	snv_base	snv_ct
#D7S820	M	AATCTGTC	GTTAGTTC	7	11	7	83789477	A	11
#D7S820	M	AATCTGTC	GTTAGTTC	7	11	7	83789520	A	11
#D7S820	M	AATCTGTC	GTTAGTTC	8	673	7	83789477	A	668
#D7S820	M	AATCTGTC	GTTAGTTC	8	673	7	83789477	C	2
#D7S820	M	AATCTGTC	GTTAGTTC	8	673	7	83789477	N	2
#D7S820	M	AATCTGTC	GTTAGTTC	8	673	7	83789477	T	1

df_summ = pd.read_csv(summ_input,sep='\t',usecols=['str_name', 'ref_chr', 'snv_pos', 'snv_base', 'motif_cnt', 'snv_ct'],
                                  dtype={'ref_chr': str, 'motif_cnt': float})
df_summ = df_summ[df_summ.ref_chr != '-']

group_cols = ['str_name', 'ref_chr', 'snv_pos', 'snv_base', 'motif_cnt']
dfsum_tots = df_summ.groupby(group_cols, sort=True).sum()

# Add index columns back into dataframe as regular columns to allow additional groupby
for i in range(5):
  dfsum_tots[group_cols[i]] = dfsum_tots.index.get_level_values(i)

df_summary = dfsum_tots.groupby(['str_name', 'ref_chr', 'snv_pos', 'snv_base']).apply(lambda row: [list(row.motif_cnt), list(row.snv_ct)])
summary_str = df_summary.to_dict()

if debug: print "\nMixture STR dictionary:\n{0}".format(summary_str)
 
#-----------------------------------------------------------------------------#
# Extract only those entries in summary_str which correspond to an entry      # 
#  in minor_str. Ie. only report minor haplotypes identified in mixture       #
#-----------------------------------------------------------------------------#
# summary_str example: 
#   ('trf605629', '4', 60021715, 'G'): [[17.5, 18.5], [196, 2]],
#   ('trf873648_trf873649', '9', 105984015, 'A'): [[25.75, 27.75], [2, 3]], ...
# 
# minor_str example:
#   ('trf605629', '4', 60021715, 'G'): ['17.5', '19.5']

print "Starting processing of minor/mixture STRs"

# Calculate total reads per SNV position in mixture
group_cols = ['str_name', 'ref_chr', 'snv_pos']
dfsnv_rds  = df_summ.groupby(group_cols, sort=True).sum()
dfsnv_rds  = dfsnv_rds.drop('motif_cnt', 1)
if debug: print dfsnv_rds

mixture_reads = dfsnv_rds.to_dict()
if debug: print mixture_reads

# Determine minor alleles present in mixture
minor_in_mix = {}
for mkey, mvalues in sorted(minor_str.items()):
  (str_name, snv_chr, snv_pos, snv_base) = mkey
  tot_snv_reads = 0
  
  if mkey in summary_str:
    for m_allele in mvalues:
      if float(m_allele) in summary_str[mkey][0]:
        idx     = summary_str[mkey][0].index(float(m_allele))
        str_rds = summary_str[mkey][1][idx]
        tot_snv_reads = mixture_reads['snv_ct'][(str_name, snv_chr, snv_pos)]

        if mkey not in minor_in_mix:
          minor_in_mix[mkey] = [[m_allele], [str_rds], tot_snv_reads]
        else:
          minor_in_mix[mkey][0].append(m_allele)
          minor_in_mix[mkey][1].append(str_rds) 

if debug:
  print "Minor_in_Mix"
  for nkey, nvalues in sorted(minor_in_mix.items()):
    print nkey, nvalues

#-----------------------------------------------------------------------------#
# Format minor haplotypes in mixture, for output                              #
#-----------------------------------------------------------------------------#
last_chr = 'N'; last_pos = 0; last_str = 'None'; last_base = 'N'
last_str_snv = 'None-N-000000'; snv_alleles = {}

haplo_csv = csv.writer(haplo_output, dialect='tab_delim')
haplo_hdgs = ['STR Name', 'SNV Chr', 'SNV Pos', 'SNV Rds', 'SNV Bases', 'STR Alleles', 'Haplotypes', 'Reads', 'Read Pcts']
haplo_csv.writerow(haplo_hdgs)

for nkey, nvalues in sorted(minor_in_mix.items()):
  # Accumulate all records for this STR (all SNV bases)
  (this_str, this_chr, this_pos, this_base) = nkey
  
  str_snv_pos = ('-').join([this_str, this_chr, str(this_pos)])
  
  if last_chr != 'N' and (this_chr != last_chr or this_pos != last_pos or this_str != last_str):
    # Have all data for previous SNV position, so determine and write out haplotype reads
    if debug: print "Generating output row for {0}, {1}, {2}".format(last_str, last_chr, last_pos)
    if debug: print snv_alleles[last_str_snv]
   
    haplo_csv.writerow(format_haplo_csv((last_str, last_chr, last_pos), snv_alleles[last_str_snv]))
 
    # Then reset STR/SNV accumulators
    snv_alleles = {}
    snv_alleles[str_snv_pos] = extract_haplotypes(this_base, nvalues)
  
  else:
    if (this_str == last_str and this_chr == last_chr and this_pos == last_pos):
      # Keep accumulating data for this STR/SNV position combination
      # Eg: snv_alleles:  {'ASTE1-3-130733039': [('C',11,889,'0.48',1825),('T',11,926,'0.51',1825)]}
      snv_alleles[str_snv_pos] += extract_haplotypes(this_base, nvalues)
    else:
      snv_alleles[str_snv_pos] = extract_haplotypes(this_base, nvalues)
      
  (last_str, last_chr, last_pos, last_base) = nkey
  last_str_snv = str_snv_pos
    
# Process last STR-SNV
haplo_csv.writerow(format_haplo_csv((last_str, last_chr, last_pos), snv_alleles[last_str_snv]))

# Close input/output files
summ_input.close()
minor_input.close()
haplo_output.close()
