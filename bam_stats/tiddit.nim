import read_bam
import vcf_header
import coverage_analysis
import extract_signals
import hts
import os
import strutils

#collect input
let
  mapq:uint8=10
  bin_size:int=1000
  output_prefix:string = "output"
  bam_file_name:string=paramStr(1)

let version = "3.0.0"

#collect library statistics
let (median_ins,ins_treshold,library_type,avg_read_length)= bam_stats(paramStr(1),50000000,mapq)

var b:Bam
open(b, bam_file_name)

#generate coverage and signal data structure
var (coverage_bins,spanning_read_bins,discordant_pair_bins,quality_bins)= make_coverage_structure(b,bin_size)
var (discordant_pairs, discordant_pairs_read_name, split_reads, split_reads_read_name)=generate_signal_structure(b.hdr)

var total_reads:uint64=0
var genome_size:uint64=0

#read the bam file
var signals:string=""

for record in b:
  if record.flag.dup:
    continue
  elif record.flag.unmapped:
    continue

  if not record.flag.supplementary and not record.flag.secondary:
    total_reads+=1
 
  #check for sv signal (split read or discordant pair
  signals.add(sv_signals(record,ins_treshold,mapq,discordant_pairs,discordant_pairs_read_name,split_reads,split_reads_read_name))
  #update coverage data
  discard update_genomic_bins(record,coverage_bins,spanning_read_bins,discordant_pair_bins,quality_bins,bin_size)

for chrom in targets(b.hdr):
  genome_size+= uint64(chrom.length)

var coverage:float

if (genome_size < 100000000):
  coverage=int(total_reads)/int(genome_size)
else:
  coverage=int(total_reads div 1000)/int(genome_size div 1000)

coverage=coverage*avg_read_length

let library_stats="##LibraryStats=TIDDIT-" & version & " Coverage=" & $coverage & " AverageReadLength=" & $avg_read_length & " MedianInsertSize=" & $median_ins & " 99.5PercentileInsertSize=" & $ins_treshold & " LibraryType=" & library_type
let output_vcf_header = print_header(library_stats,version,b.hdr)
writeFile( output_prefix & ".vcf", output_vcf_header & signals.strip())

writeFile( output_prefix & ".wig" ,print_wig(b,bin_size,coverage_bins,spanning_read_bins,discordant_pair_bins,quality_bins) )
