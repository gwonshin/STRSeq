#!/usr/bin/python

# File: str_ctlen_genotype.py
# Name: Sue Grimes
# Desc: Script summarizes STR genotypes/alleles into final output format
#
# 10/8/2015: Original version

import os, sys, csv, imp, pysam, numpy as np, pandas as pd, msi_str as msi

script_name = os.path.basename(__file__)
user_home   = os.path.expanduser("~")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if len(sys.argv) < 3:
  print "Usage: ", script_name, "<str_summary> <probe_rdcts> <str_info>"
  sys.exit(1)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Open/initialize output files and general variables                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
csv.register_dialect('tab_delim', delimiter='\t', doublequote=False, quotechar='', lineterminator='\n', escapechar='', quoting=csv.QUOTE_NONE)

summ_input  = msi.open_file(sys.argv[1], 'r')
rdct_input  = msi.open_file(sys.argv[2], 'r')
str_input   = msi.open_file(sys.argv[3], 'r')
final_fn    = sys.argv[1].replace('_summary', '_final', 1)
rpt_output  = msi.open_file(final_fn, 'w')

str_csv    = csv.DictReader(str_input, dialect='tab_delim')

FLANK_SIZE = msi.FLANK_SIZE
ALLELE2_MIN_PCT = msi.ALLELE2_MIN_PCT

print "\n**Running {0}, with summary input: {1}, probe counts: {2}".format(script_name, sys.argv[1], sys.argv[2])
print "STR file is: {0}, flank size is: {1}".format(sys.argv[3], FLANK_SIZE)
 
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
# Read STR read counts from input file                                        #
#-----------------------------------------------------------------------------#
# Probe read ct example:
# probe_nr	chromosome	probe_start_pos	str_name	strand	probe_reads	motif_reads	full_str_rds
# 1			1			3584977			trf2094		M		1248		366			111
# 2			1			3585231			trf2094		P		1872		381			204
df_cts = pd.read_csv(rdct_input,sep='\t',usecols=['str_name', 'probe_reads', 'motif_reads', 'full_str_rds'])
dfct_tots = df_cts.groupby(['str_name'], sort=True).sum()
str_probe_cts = dfct_tots.to_dict(orient='dict')

#-------------------------------------------------------------------------------#
# Determine genotype(s) for each STR, by assessing read counts per motif repeat # 
#   in summary input file                                                       #
#-------------------------------------------------------------------------------#
# Summary input example:
# STRName Strand	Motif	Min Rpts	Probe Rds	5pr Flank	3pr Flank		STR Len	Motif#	STR Rds
#  AR		M		GA		3			4818	GCCTCAATGAACTGG	CAGCTTGTACACGTG	8		4		1964
#  AR		P		GA		3			4199	GCCTCAATGAACTGG	CAGCTTGTACACGTG	8		4		2208
df_summ    = pd.read_csv(summ_input,sep='\t',usecols=['STR Name', 'STR Len', 'Motif#', 'STR Rds'])
df_summ.columns = ['str_name', 'str_len', 'motif_ct', 'str_rds']

group_cols = ['str_name', 'str_len', 'motif_ct']
dfsum_tots = df_summ.groupby(group_cols, sort=True).sum()

# Add index columns back into dataframe as regular columns to allow additional groupby
for i in range(3):
  dfsum_tots[group_cols[i]] = dfsum_tots.index.get_level_values(i)

df_summary = dfsum_tots.groupby('str_name').apply(lambda row: [list(row.motif_ct), list(row.str_rds)])
summary_str = df_summary.to_dict()

#-----------------------------------------------------------------------------#
# Determine alleles, and write final output                                   #
#-----------------------------------------------------------------------------#
final_csv = csv.writer(rpt_output, dialect='tab_delim')
final_csv.writerow(['STR Name', 'Motif', 'Min Rpts', 'Probe Rds', '5pr Flank', '3pr Flank', 'STR TotRds', 'Motif Ct', 'STR Rds', 'Allele(s)'])

for str_name, svalues in sorted(summary_str.items()):
  str_nm_info  = str_info[str_name]
  
  probe_rds = str_probe_cts['probe_reads'][str_name] if str_name in str_probe_cts['probe_reads'] else 0
  motif_rds   = str_probe_cts['motif_reads'][str_name] if str_name in str_probe_cts['motif_reads'] else 0
  full_str_rds = str_probe_cts['full_str_rds'][str_name] if str_name in str_probe_cts['full_str_rds'] else 0
  
  motif_cts = ', '.join("{0:.5g}".format(motif_rpts) for motif_rpts in svalues[0])
  str_reads = ', '.join(str(str_rds) for str_rds in svalues[1])
  [status, str_alleles] = msi.determine_alleles(svalues[1], svalues[0])
  allele_string = ', '.join("{0:.5g}".format(float(str_allele)) for str_allele in str_alleles)
  
  final_csv.writerow([str_name, str_nm_info['str_motif'], str_nm_info['min_repeats'], probe_rds, str_nm_info['flanking_5pr'],  
                       str_nm_info['flanking_3pr'], full_str_rds, motif_cts, str_reads, allele_string])

# Close output files
rpt_output.close()
