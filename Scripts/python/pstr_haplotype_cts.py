
#!/usr/bin/python

# File: pstr_haplotype_cts.py
# Name: Sue Grimes
# Desc: Script takes final STR/SNV information and reformats for haplotype read ct normalization
#
# 12/17/2015: Original version
# 1/8/2016:   Allow specifying reporting of all alleles, or just major alleles
# 2/11/2016:  Modify to use str allele thresholds when determining haplotypes
#             Use pandas for some of the data structure manipulation/sorting
# 2/16/2016:  Modify to handle case where same read counts for top 2 alleles,
#             Sort by str repeat count (ascending) to pick highest repeat count as major.

import os, sys, csv, imp, pandas as pd, msi_str as msi
from decimal import Decimal
from collections import OrderedDict

script_name = os.path.basename(__file__)
user_home   = os.path.expanduser("~")

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if len(sys.argv) < 2:
  print "Usage: ", script_name, "<str_snv_final_file> [major|all]"
  sys.exit(1)

alleles_to_report = 'all'
if len(sys.argv) > 2 and sys.argv[2] == 'major':
  alleles_to_report = 'major'

debug = False
if len(sys.argv) > 3 and sys.argv[3] == 'debug':
  debug = True 

MINOR_BASE_MIN = 0.15

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Open/initialize output files and general variables                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
csv.register_dialect('tab_delim', delimiter='\t', doublequote=False, quotechar='', lineterminator='\n', quoting=csv.QUOTE_NONE)

strsnv_input  = msi.open_file(sys.argv[1], 'r')
strsnv_csv = csv.DictReader(strsnv_input, dialect='tab_delim') 

haplo_fbase = sys.argv[1].split('/')[-1].split('.')[0]
haplo_output = msi.open_file(haplo_fbase + '.haplotype_cts_' + alleles_to_report + '.txt', 'w')
haplo_csv = csv.writer(haplo_output, dialect='tab_delim')

print "\n**Running {0}, with STR-SNV input: {1}".format(script_name, sys.argv[1])

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# General methods                                                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
def extract_haplotypes(crow):
  # Returned tuples: [('C','11',889), ('C','10',96), ('C','9',9), ('C','12',5)]
  snv_base      = crow['SNV Allele']
  str_alleles   = [s_allele.lstrip() for s_allele in crow['Motif Rpts'].split(',')]
  snv_reads     = [int(rds) for rds in crow['Motif Rpt Rds'].split(',')]
  
  haplo_tuples = []
  for i, allele in enumerate(str_alleles):
    haplo_tuples.append((snv_base, allele, snv_reads[i]))
  if debug: print "Haplotype tuples (all): {0}".format(haplo_tuples)
  return haplo_tuples
  
def determine_haplotypes(snv_allele_rds):
# snv_allele_rds: [(C,11,889),(T,11,926),(T,10,454)]
  df = pd.DataFrame(snv_allele_rds, columns=['base', 'str_rpt', 'rds'])
  df_bases = df.groupby(['base']).sum().reset_index()
  df_bases.sort_values(by='rds', ascending=False, inplace=True)
  major_bases = major_snv_bases(df_bases)
  
  df_repeats = df.groupby(['str_rpt']).sum().reset_index()
  df_repeats.sort_values(by='rds', ascending=False, inplace=True) 
  all_rpt_alleles =  df_repeats.values.tolist()
  
  status = 'None'
  if len(all_rpt_alleles) > 1:
    for potential_minor_allele in all_rpt_alleles[1:]:
      [status, major_alleles] = msi.check_threshold(all_rpt_alleles[0], potential_minor_allele,1) 
      if status == 'Done':  break    
  
  if status != 'Done':
    major_alleles = [all_rpt_alleles[0][0]]
  
  nr_bases   = len(major_bases)
  nr_alleles = len(major_alleles)
  
  if nr_alleles == 1 and nr_bases == 1:
    return [htuple for htuple in snv_allele_rds if htuple[0] == major_bases[0] and htuple[1] == major_alleles[0]]
  elif nr_alleles == 1 and nr_bases == 2:
    return [htuple for htuple in snv_allele_rds if htuple[0] in major_bases and htuple[1] == major_alleles[0]]
  elif nr_alleles == 2 and nr_bases == 1:
    return [htuple for htuple in snv_allele_rds if htuple[0] == major_bases[0] and htuple[1] in major_alleles]
  elif nr_alleles == 2 and nr_bases == 2:
    potential_alleles = [htuple for htuple in snv_allele_rds if htuple[0] in major_bases and htuple[1] in major_alleles]
    return major_haplotypes(potential_alleles)
  else:
    return snv_allele_rds
  
def major_snv_bases(df_bases):
  # Need to add logic for case where >2 bases, and 2nd & 3rd base #reads are equal?
  # Eg. [['A',10], ['G',3], ['T',3]]
  bases_reads = df_bases.values.tolist()
  if len(bases_reads) == 1 or (bases_reads[1][1] < (bases_reads[0][1] * MINOR_BASE_MIN)):
    return [bases_reads[0][0]]
  else:
    return [bases_reads[0][0], bases_reads[1][0]]
 
def major_haplotypes(potential_alleles):
  dfa = pd.DataFrame(potential_alleles, columns=['base', 'str_rpt', 'rds'])
  dfg = dfa.groupby(['base','str_rpt']).sum().reset_index()
  dfg.sort_values(by=['rds', 'str_rpt'], ascending=False, inplace=True)
  alleles_srted = dfg.values.tolist()
  final_alleles = [alleles_srted[0]]
  
  for i in range(1,len(alleles_srted)):
    if alleles_srted[i][0] != final_alleles[0][0] and alleles_srted[i][1] != final_alleles[0][1]:
      final_alleles += [alleles_srted[i]]
      break
  return final_alleles  

def format_haplo_csv(snv_allele_rds):
# snv_allele_rds: [(C,11,889),(T,11,926),(T,10,454)]
  alleles_rds_srted = sorted(snv_allele_rds, key=lambda x: int(x[2]), reverse=True)
  haplotypes  = [('-').join([base, str_allele]) for (base, str_allele, allele_rds) in alleles_rds_srted]
  
  snv_bases   = [base for (base, str_allele, allele_rds) in alleles_rds_srted]
  snv_bases_uniq = list(OrderedDict.fromkeys(snv_bases))
  
  str_alleles = [str_allele for (base, str_allele, allele_rds) in alleles_rds_srted]
  str_alleles_uniq = list(OrderedDict.fromkeys(str_alleles))
  
  haplo_rds   = [allele_rds for (base, str_allele, allele_rds) in alleles_rds_srted]
  tot_snv_reads = sum(haplo_rds)
  hap_rd_pcts  = [(allele_rds * 100.00/tot_snv_reads) for (base, str_allele, allele_rds) in alleles_rds_srted]
  
  return [last_str, last_chr, last_pos, tot_snv_reads, ','.join(snv_bases_uniq), ','.join(str_alleles_uniq),
            ','.join(haplotypes), ','.join([str(rds) for rds in haplo_rds]),
            ','.join(["{0:.2f}".format(hap_rd_pct) for hap_rd_pct in hap_rd_pcts])] 
  
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Main program logic                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# STR input example (some unused columns omitted for brevity/clarity):
# STR name	Motif	Ref Chr	SNV Pos		SNV Allele	SNV Reads	Motif Rpts	Motif Rpt Rds	STR Allele(s)
# ATM_2		AT		11		108143311		A		32			4			32				4
# BRCA1_1	CA		17		41215825		C		1181		3, 4		1, 1180			4
# BRCA1_1	CA		17		41215825		T		1130		4			1130			4
# BRCA1_2	T		17		41245466		A		1452		7, 8, 9		5, 1445, 2		8
# BRCA1_2	T		17		41245466		G		1298		7, 8, 10	3, 1294, 1		8

# Output (major alleles) example:
# STR Name	SNV Chr	SNV Pos	SNV Reads	SNV Bases	STR Alleles	Haplotypes	Reads
# ATM_2		11		108143311	32			A			4		A-4			32
# BRCA1_1	17		41215825	2311		C,T			4		C-4,T-4		1180,1130
# BRCA1_2	17		41245466	2750		A,G			8		A-8,G-8		1445,1294

last_chr = 'N'; last_pos = 0; last_str = 'None'; 
snv_alleles = {}
haplo_csv.writerow(['STR Name', 'SNV Chr', 'SNV Pos', 'SNV Reads', 'SNV Bases', 'STR Alleles', 'Haplotypes', 'Reads', 'Read Pcts'])

for crow in strsnv_csv:
  # Accumulate all records for this STR (all SNV bases)
  this_chr = crow['Ref Chr']
  this_pos = crow['SNV Pos']
  this_str = crow['STR name']
  str_snv_pos = ('-').join([this_str, this_chr, str(this_pos)])
  
  if last_chr != 'N' and (this_chr != last_chr or this_pos != last_pos or this_str != last_str):
    # Have all data for previous SNV position, so determine and write out haplotype reads
    print "Generating output row for {0}, {1}, {2}".format(last_str, last_chr, last_pos)
    if debug: print snv_alleles[last_str_snv]
   
    final_alleles = snv_alleles[last_str_snv]
    if alleles_to_report == 'major':
      final_alleles = determine_haplotypes(snv_alleles[last_str_snv])
    if final_alleles and len(final_alleles) > 0:
      haplo_csv.writerow(format_haplo_csv(final_alleles))
    else:
      print "*Error - empty array.  SNV alleles: {0}".format(snv_alleles[last_str_snv])
 
    # Then reset STR/SNV accumulators
    snv_alleles = {}
    snv_alleles[str_snv_pos] = extract_haplotypes(crow)
  
  else:
    if (this_str == last_str and this_chr == last_chr and this_pos == last_pos):
      # Keep accumulating data for this STR/SNV position combination
      # Eg: snv_alleles:  {'ASTE1-3-130733039': [(C,11,889),(T,11,926),(T,10,454)]}

      snv_alleles[str_snv_pos] += extract_haplotypes(crow)
    else:
      snv_alleles[str_snv_pos] = extract_haplotypes(crow)
      
  last_chr = this_chr
  last_pos = this_pos
  last_str = this_str
  last_str_snv = str_snv_pos
    
# Process last STR-SNV
final_alleles = snv_alleles[last_str_snv]
if alleles_to_report == 'major':
  final_alleles = determine_haplotypes(snv_alleles[last_str_snv])
if len(final_alleles) > 0:
  haplo_csv.writerow(format_haplo_csv(final_alleles))
else:
  print "*Error - empty array.  SNV alleles: {0}".format(snv_alleles[last_str_snv])

# Close input/output files
strsnv_input.close()
haplo_output.close()
