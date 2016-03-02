#!/usr/bin/python

# File: pstr_merge_str_snv.py
# Name: Sue Grimes
# Desc: Script merges STR information from R1, and SNV information from R2,
#          based on query name
#
# 6/?/2015:  Original version
# 10/6/2015: Adding heading comments

import os, sys, csv, imp, pandas as pd, msi_str as msi

script_name = os.path.basename(__file__)
user_home   = os.path.expanduser("~")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid command line arguments                                      #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if len(sys.argv) < 2:
  print "Usage: ", script_name, " <R1_STR_file> <R2_SNV_file>"
  sys.exit(1)

str_ifile = sys.argv[1]
if not os.path.isfile(str_ifile) or not os.access(str_ifile, os.R_OK):
  print "Unable to open STR file for input: ", str_ifile
  sys.exit(1)

snv_ifile = sys.argv[2]
if not os.path.isfile(snv_ifile) or not os.access(snv_ifile, os.R_OK):
  print "Unable to open SNV file for input: ", snv_ifile
  sys.exit(1)
  
fnbase = str_ifile.split('/')[-1].split('.')[0]

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Open/initialize output files and general variables                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
csv.register_dialect('tab_delim', delimiter='\t', doublequote=False, quotechar='', lineterminator='\n', escapechar='', quoting=csv.QUOTE_NONE)
dtl_ofile = fnbase + '.STR_SNV.detail.txt'
summ_ofile = fnbase + '.STR_SNV.summary.txt'

print "\n**Running {0}, with R1 STR: {1}, R2 SNV:{2}".format(script_name, str_ifile, snv_ifile)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Main program logic                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#-----------------------------------------------------------------------------#
# Merge STR and SNV files on query name                                       #
#-----------------------------------------------------------------------------#
df_str = pd.read_csv(str_ifile,sep='\t',dtype={'chrom': str})
df_snv = pd.read_csv(snv_ifile,sep='\t',dtype={'ref_chr': str})

df = pd.merge(df_str,df_snv,on=('qname','probe_nr'),how='left')
# Since not all STR reads will have an SNV in the matching R2, and we are doing a left join,
#   integer columns which have missing values, will be converted to floating point, use
#   format command to remove .0
# Use fillna on the snv_base column since that will be a key field and 'na' columns are ignored
#   in groupby
df['ref_chr'] = df['ref_chr'].fillna('-')
df['snv_pos'] = df['snv_pos'].fillna(0)
df['snv_pos'] = df['snv_pos'].map('{:.0f}'.format)
df['snv_base'] = df['snv_base'].fillna('-')

# Summarize, counting reads by STR repeat allele
key_cols=['probe_nr','str_name','strand', 'bases_5pr','bases_3pr','motif_cnt']
dfstr_cts = pd.DataFrame({'str_ct': df_str.groupby(key_cols, sort=True)['seq'].count()}).reset_index()
print dfstr_cts.describe
#dfsumm['qname'].to_csv(summ_ofile,index=True,sep='\t',header=['read_ct'])

# Summarize merged file, counting SNV alleles by STR repeat allele
key_cols=['probe_nr','str_name','strand', 'bases_5pr','bases_3pr','motif_cnt', 'ref_chr', 'snv_pos', 'snv_base']
dfsnv_cts = pd.DataFrame({'snv_ct': df.groupby(key_cols, sort=True)['seq'].count()}).reset_index()
print dfsnv_cts.describe

# Merge STR and SNV allele counts, and summarize by probe#
dfsumm = pd.merge(dfstr_cts,dfsnv_cts,on=('probe_nr','str_name','strand','bases_5pr','bases_3pr','motif_cnt'),how='left')
print dfsumm
summ_cols = ['probe_nr', 'str_name', 'strand', 'bases_5pr', 'bases_3pr', 'motif_cnt', 'str_ct', 'ref_chr', 'snv_pos', 'snv_base', 'snv_ct']
dfsumm[summ_cols].to_csv(dtl_ofile,index=False,sep='\t')

# Summarize by STR name
summ_cols = ['str_name', 'strand', 'bases_5pr', 'bases_3pr', 'motif_cnt', 'str_ct', 'ref_chr', 'snv_pos', 'snv_base', 'snv_ct']
dfsumm[summ_cols].to_csv(summ_ofile,index=False,sep='\t')

