#!/usr/bin/python

# File: str_flank_alleles.py
# Name: Sue Grimes
# Desc: Script outputs ref and alt alleles in STR flanking regions
#
# 10/6/2015: Original version

import os, sys, csv, imp, msi_str as msi

script_name = os.path.basename(__file__)
user_home   = os.path.expanduser("~")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid command line arguments                                      #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if len(sys.argv) < 2:
  print "Usage: ", script_name, "<flank_variants_file>"
  sys.exit(1)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Open/initialize output files and general variables                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
csv.register_dialect('tab_delim', delimiter='\t', doublequote=False, quotechar='', lineterminator='\n', escapechar='', quoting=csv.QUOTE_NONE)
snv_outfn  = sys.argv[1].split(".")[0] + '.flank_alleles.txt'

in_csv = csv.reader(msi.open_file(sys.argv[1],'r'), dialect='tab_delim')
out_csv = csv.writer(msi.open_file(snv_outfn, 'w'), dialect='tab_delim')

print "\n**Running {0}".format(script_name)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Main program logic                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#-----------------------------------------------------------------------------#
# Read STR flanking SNV information and write alternate flanking alleles      #
#-----------------------------------------------------------------------------#
out_csv.writerow(['Chr', 'SNVPos', 'Ref', 'Alt', 'TYPE', 'GT', 'AF', 'STRName', '5or3pr', 'FlankStart', 'FlankEnd'])
for srow in in_csv:
  gt_parsed = srow[9].split(':')
  info_parsed = srow[7].split(';')
  var_type = [fld[5:] for fld in info_parsed if fld[0:5] == 'TYPE=']
  
  if var_type[0] in ['snp', 'snp,snp']:
    genotype = gt_parsed[0]
    allele_freq = [fld[3:] for fld in info_parsed if fld[0:3] == 'AF=']
    out_csv.writerow([srow[0], srow[1], srow[3], srow[4], var_type[0], genotype, allele_freq[0], srow[14], srow[10], srow[12], srow[13]])

