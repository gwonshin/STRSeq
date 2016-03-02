#!/usr/bin/python

# File: pstr_genotyping.py
# Name: Sue Grimes
# Desc: Script takes summary STR/SNV information and reformats into final output
#
# 6/12/2015: Original version
# 10/6/2015: Use msi.open_file method for input/output files
# 10/19/2015: Use pandas to create summary file
# **NOTE: Some input file and dictionary examples need to be modified **
# 11/18/2015: Use new method (from GiWon) for determination of alleles

import os, sys, csv, imp, pysam, numpy as np, pandas as pd, msi_str as msi
from decimal import Decimal

script_name = os.path.basename(__file__)
user_home   = os.path.expanduser("~")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if len(sys.argv) < 4:
  print "Usage: ", script_name, "<str_snv_summary_file> <str_info> <probe_cts> [debug]"
  sys.exit(1)

debug = False
if len(sys.argv) > 4 and sys.argv[4] == 'debug':
  debug = True 

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Open/initialize output files and general variables                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
csv.register_dialect('tab_delim', delimiter='\t', doublequote=False, quotechar='', lineterminator='\n', quoting=csv.QUOTE_NONE)

summ_input  = msi.open_file(sys.argv[1], 'r')
str_input   = msi.open_file(sys.argv[2], 'r')
probe_input = msi.open_file(sys.argv[3], 'r')

str_csv = csv.DictReader(str_input, dialect='tab_delim')
summ_fbase = sys.argv[1].split('/')[-1].split('.')[0] 

final_output = msi.open_file(summ_fbase + '.STR_SNV.final.txt', 'w')
final_csv = csv.writer(final_output, dialect='tab_delim')

FLANK_SIZE = msi.FLANK_SIZE
ALLELE2_MIN_PCT = msi.ALLELE2_MIN_PCT

print "\n**Running {0}, with STR/SNV input: {1}".format(script_name, sys.argv[1])
print "Parameters: Flank size: {0}".format(FLANK_SIZE)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# General methods                                                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#


#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Main program logic                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#-----------------------------------------------------------------------------#
# Create str_info dictionary, from input csv file                             #
#-----------------------------------------------------------------------------#
# STR input example:
#  STRName	Motif	MinRepeat	5PrEnd	5PrFlank		3PrStart	3PrFlank
#  PIK3CA	A		5			6273772 GCTGAAGCAATCAGG	6273799		CTCGAAGTATGTTGC
str_info   = msi.dict_from_csv(str_csv, 'str_info')

#-----------------------------------------------------------------------------#
# Read probe/str counts file, and load into dictionary                        #
#-----------------------------------------------------------------------------# 
#probe_nr	chromosome	probe_start_pos	str_name	strand	probe_reads	full_str_rds
#664	3	189420981	trf583174	M	419	14
#665	3	189421335	trf583174	P	2760	697
#417	2	152141668	trf451304	M	3365	546
#418	2	152141823	trf451304	M	1269	554
#
df_cts = pd.read_csv(probe_input,sep='\t',usecols=['str_name', 'probe_reads', 'motif_reads', 'full_str_rds'])
dfct_tots = df_cts.groupby(['str_name'], sort=True).sum()
str_probe_cts = dfct_tots.to_dict(orient='dict')
if debug: print "\nSTR Probe Counts:\n{0}".format(str_probe_cts)

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
if debug: print "\nSummary STR dictionary:\n{0}".format(summary_str)
 
#-----------------------------------------------------------------------------#
# Reformat into one line per STRname/SNV with STR alleles as columns          #
#-----------------------------------------------------------------------------#
# summary_str example: 
#   ('trf605629', '4', 60021715, 'G'): [[17.5, 18.5], [196, 2]],
#   ('trf873648_trf873649', '9', 105984015, 'A'): [[25.75, 27.75], [2, 3]], ...
# 
final_csv = csv.writer(final_output, dialect='tab_delim')
final_hdgs = ['STR name', 'Motif', 'Min Rpts', 'Probe Rds', 'Rds w/motif', 'Flank-5pr', 'Flank-3pr', 'Full STR Rds', 'Ref Chr', 'SNV Pos',
               'SNV Allele', 'SNV Reads', 'Motif Rpts', 'Motif Rpt Rds', 'STR Allele(s)']
final_csv.writerow(final_hdgs)

for summ_key, svalues in sorted(summary_str.items()):
  print summ_key, svalues
  (str_name, snv_chr, snv_pos, snv_base) = summ_key  
  str_nm_info = str_info[str_name]
 
  probe_rds = str_probe_cts['probe_reads'][str_name] if str_name in str_probe_cts['probe_reads'] else 0
  motif_rds   = str_probe_cts['motif_reads'][str_name] if str_name in str_probe_cts['motif_reads'] else 0
  full_str_rds = str_probe_cts['full_str_rds'][str_name] if str_name in str_probe_cts['full_str_rds'] else 0

  str_reads = ', '.join(str(str_rds) for str_rds in svalues[1])
  motif_rpts = ', '.join("{0:.5g}".format(motif_cts) for motif_cts in svalues[0])
  [status, str_alleles] = msi.determine_alleles(svalues[1], svalues[0])
  allele_string = ', '.join("{0:.5g}".format(float(str_allele)) for str_allele in str_alleles)
  
  final_csv.writerow([str_name, str_nm_info['str_motif'], str_nm_info['min_repeats'], probe_rds, motif_rds, 
                       str_nm_info['flanking_5pr'], str_nm_info['flanking_3pr'], full_str_rds,
                       snv_chr, snv_pos, snv_base, sum(svalues[1]), motif_rpts, str_reads, allele_string])

# Close input/output files
summ_input.close()
str_input.close()
final_output.close()
