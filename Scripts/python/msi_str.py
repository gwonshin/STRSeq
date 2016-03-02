#!/usr/bin/env python
import sys, os, subprocess, re, string, itertools
import numpy as np

FLANK_SIZE = int(os.getenv('FLANK_SIZE', 15))
MAX_FLANK_MISMATCH = 2
ALLELE2_MIN_PCT = 0.5
THRESHOLD_VALS = [0.45, 0.35, 0.15, 0.02]

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define reusable modules to be imported/called from STR scripts              #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
def open_file(fname, rw):
  rw_type = {'r': 'input', 'w': 'output'}
  try:
    return open(fname, rw)
  except:
    input_output = 'unknown' if rw not in rw_type else rw_type[rw]
    print "Unable to open {0} file: {1}".format(input_output, fname)
    sys.exit(1)

def dict_from_csv(csv_in, ftype):
  cdict = {}
  if ftype == 'str_info':
    for crow in csv_in:
      str_name = crow['STRName']
      cdict[str_name] = {'str_motif': crow['Motif'], 'min_repeats': int(crow['MinRepeat']),
                         'flanking_5pr': flank_seq(crow['5PrFlank'],'5pr'), '5r_end_pos': crow['5PrEnd'],
                         'flanking_3pr': flank_seq(crow['3PrFlank'],'3pr'), '3pr_start_pos': crow['3PrStart']}
  
  elif ftype == 'probe_info':
    for crow in csv_in:
      probe_nr = int(crow['Probe'])
      cdict[probe_nr] = {'chromosome': crow['Chr'], 'probe_start_pos': int(crow['StartPos']), 
                         'str_name': crow['STRName'], 'strand': crow['Strand'], 'probe_reads': 0,
                         'motif_reads': 0, 'full_str_rds': 0}

  elif ftype == 'probe_cts':
    for crow in csv_in:
      str_name_strand = crow['str_name'] + '-' + crow['strand']
      cdict[str_name_strand] = {'probe_rds': int(crow['probe_reads']), 'motif_reads': int(crow['motif_reads'])}

  else:
    raise ValueError('Invalid file type for dict_from_csv method', ftype)
  return cdict

def increment_dict_ct(dict, dkey):
  if dkey in dict.keys():
    dict[dkey] += 1
  else:
    dict[dkey] = 1

def flank_seq(flankseq, pr5or3):
  if pr5or3 == '5pr':
    return flankseq[0-FLANK_SIZE:]
  else:
    return flankseq[0:FLANK_SIZE]

def motif_search(seq, rlen, motif, min_repeat):
  multi_motif  = ''.join([motif]*min_repeat)
  motif_length = len(motif)
  str_values   = {'motif_cnt': 0}
  
  if seq.find(multi_motif) >= 0:
    # For a motif 'GATA' and min repeat of 3, regex string will be: r'((GATA){3,})'
    regex_string = r'((%s){%d,})' % (motif, min_repeat)
    motif_match = re.search(regex_string, seq)
    matched_string = motif_match.group(1)
    motif_idx   = seq.find(matched_string) #Get starting position of matched string
    motif_cnt   = len(matched_string) / motif_length #Calculate number of repeats for motif

    # Check for possible truncation of STR at beginning(5') or end(3') of read
    flanking_5pr_start = max(0, motif_idx - FLANK_SIZE)
    flanking_5pr_end   = max(0, motif_idx)
    len_flanking_5pr   = flanking_5pr_end - flanking_5pr_start

    flanking_3pr_start = min(rlen, motif_idx + len(matched_string))
    flanking_3pr_end   = min(rlen, motif_idx + len(matched_string) + FLANK_SIZE)
    len_flanking_3pr   = flanking_3pr_end - flanking_3pr_start
      
    trunc_flag = ''
    if len_flanking_5pr >= FLANK_SIZE and len_flanking_3pr >= FLANK_SIZE:
      trunc_flag = 'ok'
    else:
      if len_flanking_5pr < FLANK_SIZE:
        trunc_flag = '*b'
      if len_flanking_3pr < FLANK_SIZE:
        trunc_flag = trunc_flag + '*e'
      
    # Return output dictionary (hash)
    # All values are relative to 0 where 0 is beginning of genomic sequence
    str_values = {'motif_idx': motif_idx, 'motif_cnt': motif_cnt, 'matched_string': matched_string, 'f5pr_beg': flanking_5pr_start, 'f5pr_end': flanking_5pr_end,
                  'f3pr_beg': flanking_3pr_start, 'f3pr_end': flanking_3pr_end, 'trunc_flag': trunc_flag}
  return str_values 

def rev_complement(seq):
  base_complement = string.maketrans('ACTGN.', 'TGACNN')
  return seq.translate(base_complement)[::-1]

def chr_pos_motif(align_start, motif_idx, strand):
  motif_pos = align_start + motif_idx 
  return motif_idx if align_start == -1 else motif_pos
      
def snv_ref_alt(ref_base, alt_alleles, allele_freq):
  # Determine alternate alleles and their frequencies
  # Return ref base/frequency, along with all alternate alleles and their frequencies
  # Format of returned array: [('G',0),('T',0.5),('A',0.5)] if heterozygous for two alt alleles
  alt_bases = alt_alleles.split(',')
  alt_freq  = allele_freq.split(',')
  if len(ref_base) == 1:
    ref_alt_bases = [ref_base] + alt_bases
  else:
    ref_alt_bases = snv_bases([ref_base] + alt_bases)
  af_ratios = [float(freq) for freq in alt_freq]
  all_freq  = [1 - sum(af_ratios)] + af_ratios
  return zip(ref_alt_bases, all_freq)
  
def snv_bases(ref_alt_bases):
  # Handle the situation where FreeBayes reports >1 base string for SNV (eg. Ref: ATATCG  Alt: ATATAG)
  seq_len = len(ref_alt_bases[0])
  # snv_pos is array of positions in string where reference is different to variant;
  # This is an SNV, so should only be one position different => just return ref/alt bases bases at first position 
  #   where ref and alt are different (ie snv_pos[0])
  snv_pos = [i for i in xrange(seq_len) if ref_alt_bases[0][i] != ref_alt_bases[1][i]]
  return [ref_alt_base[snv_pos[0]] for ref_alt_base in ref_alt_bases]

def determine_alleles(str_rds, motif_rpts):
  nr_motif_cts = len(motif_rpts)
  rd_ct_instances = dict((rd,str_rds.count(rd)) for rd in set(str_rds))
  # Example rd_ct_instances: {357: 1, 210: 1, 3: 2, 52: 1, ..}
  
  merged_arr = [[str_rds[i], motif_rpts[i]] for i in range(0, len(str_rds))]
  sorted_arr = sorted(merged_arr, reverse=True) 
  # Example sorted_arr:  [[357,19], [210,12], [52,18], [23,13], ..]
  
  i = 0; result= ['Next',[0]]
  
  #Evaluate major allele
  [rds, motif_rpt] = sorted_arr[0]
  if nr_motif_cts == 1:
    result = ['Done',[motif_rpt]]
  elif rd_ct_instances[rds] == 2:
    result = ['Done',[motif_rpt, sorted_arr[i+1][1]]]
  elif rd_ct_instances[rds] > 2:
    result = ['Done',[-1]]
  else:  # Only one major allele; Evaluate minor alleles
    while (result[0] == 'Next'):
      i += 1
      if i == nr_motif_cts:
        result = ['Done',[sorted_arr[0][1]]]
      else:
        result = evaluate_minor_alleles(rd_ct_instances, sorted_arr, i)
  return result
  
def evaluate_minor_alleles(rd_ct_instances, sorted_arr, i):
  max_rd_values    = sorted_arr[0]
  min_rd_values    = sorted_arr[-1]
  max_motif_rpt    = max_rd_values[1]
  min_reads        = min_rd_values[0]
  shoulder_repeats = [max_motif_rpt-1, max_motif_rpt+1]
  [rds, motif_rpt] = sorted_arr[i]
  result = ['Next',[0]]
  
  if abs(motif_rpt - max_motif_rpt) <= 1.5 or local_max(sorted_arr,i) == True:
    if rd_ct_instances[rds] == 1:
      result = check_threshold(max_rd_values, sorted_arr[i])
    else:
      potential_minor_alleles = [[mrds, mrpt_ct] for [mrds, mrpt_ct] in sorted_arr if mrds == rds]
      best_minor = get_best_minor(potential_minor_alleles, max_rd_values)
      result = check_threshold(max_rd_values, best_minor[0])
      if len(best_minor) > 1 and result[0] == 'Done':
        print "Best minor", best_minor, result
        result = ['Done',[max_motif_rpt, -2]]
  return result
  
def get_best_minor(potential_minor_alleles, max_rd_values):
  max_rpt_ct         = max_rd_values[1]
  rds                = potential_minor_alleles[0][0]
  rpt_diff_from_max  = sorted([p_cts-max_rpt_ct for [p_rds, p_cts] in potential_minor_alleles])
  if max(rpt_diff_from_max) > 0:
    # Return only repeat counts to right of max which have lower threshold to pass than left of max
    best_alleles = [[p_rds, p_cts] for [p_rds, p_cts] in potential_minor_alleles if p_cts > max_rpt_ct]
  elif max(rpt_diff_from_max) == -1:
    # Keep just left shoulder which has lower threshold to pass than further left
    best_alleles = [[rds, max_rpt_ct-1]]
  else:
    # Keep everything
    best_alleles = potential_minor_alleles
  return best_alleles
  
def local_max(sorted_arr, i):
  i_rpts   = sorted_arr[i][1]
  rds_in_window = [rds for [rds, rpt_ct] in sorted_arr if rpt_ct <= i_rpts+1 and rpt_ct >= i_rpts-1]
  return True if sorted_arr[i][0] == max(rds_in_window) else False

def local_max_window2(sorted_arr, i):
  max_rpts = sorted_arr[0][1]
  i_rpts   = sorted_arr[i][1]
  if i_rpts > max_rpts:
    w_min = max(max_rpts+1, i_rpts-2)
    w_max = i_rpts+2
  else:
    w_min = i_rpts-2
    w_max = min(max_rpts-1, i_rpts+2)

  rds_in_window = [rds for [rds, rpt_ct] in sorted_arr if rpt_ct <= w_max and rpt_ct >= w_min]
  return True if sorted_arr[i][0] == max(rds_in_window) else False
  
def check_threshold(max_rd_values, i_values, idx_rds=0):
  print "Checking threshold: Max allele {0}, Test allele {1}, Reads in position: {2}".format(max_rd_values, i_values, idx_rds)
  idx_mrpt       = 1 if idx_rds == 0 else 0; 
  motif_rpt_diff = float(i_values[idx_mrpt]) - float(max_rd_values[idx_mrpt])
  motif_reads    = i_values[idx_rds]
  max_reads      = max_rd_values[idx_rds]
  
  if motif_rpt_diff < -1.0001:
    threshold_reads = max_reads * THRESHOLD_VALS[0] 
  elif motif_rpt_diff < 0:
    threshold_reads = max_reads * THRESHOLD_VALS[1]
  elif motif_rpt_diff > 1.0001:
    threshold_reads = max_reads * THRESHOLD_VALS[3] 
  else: # motif_rpt_diff > 0
    threshold_reads = max_reads * THRESHOLD_VALS[2]
  
  #print max_rd_values, i_values, motif_rpt_diff, threshold_reads
  if motif_reads > threshold_reads:
    return ['Done',[max_rd_values[idx_mrpt], i_values[idx_mrpt]]]
  else:
    return ['Next',[0]]