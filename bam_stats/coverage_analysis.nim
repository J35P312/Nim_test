import hts
import os
import sequtils
import tables
import math
import system

proc make_coverage_structure(b: Bam, bin_size: int) :Table[string,seq[float]]=
  var coverage_bins = initTable[string, seq[float]]()

  var header= b.hdr
  for line in targets(header):
    coverage_bins[line.name]= newSeqWith( int(ceil( float(line.length)/float(bin_size))) , float(0))

  return(coverage_bins)

proc update_genomic_bins(alignment: Record,genome_bins:var Table[string,seq[float]],bin_size:int): int=

 let
   first_bin=int(floor( int(alignment.start) / bin_size ))
   last_bin=int(floor(int(alignment.stop)/bin_size))
 if first_bin == last_bin:
   genome_bins[alignment.chrom][first_bin] += (1+int(alignment.stop-alignment.start))/bin_size
 else:
   for i in first_bin..last_bin:
     #echo alignment.start , " " , alignment.stop , " " , first_bin , " " , last_bin
     if i == first_bin:
       genome_bins[alignment.chrom][i] += (bin_size-(int(alignment.start)-first_bin*bin_size)+1)/bin_size
     elif i == last_bin:
       genome_bins[alignment.chrom][i] += (int(alignment.stop)-last_bin*bin_size+1)/bin_size
     else:
       genome_bins[alignment.chrom][i] += float(1)

proc print_bed(b:Bam,bin_size:int, genome_bins: Table[string,seq[float]]): int =
  for chrom in targets(b.hdr):
    for i in 0 .. genome_bins[chrom.name].len-1:
      echo chrom.name , "\t" , 1+i*bin_size, "\t" , (i+1)*bin_size ,"\t", genome_bins[chrom.name][i]

var file_name:string=paramStr(1)

var b:Bam
open(b, file_name)

let
  mapq:uint8=10
  bin_size:int=1000

var coverage_bins= make_coverage_structure(b,bin_size)

for record in b:
    if record.flag.dup:
       continue
    elif record.flag.unmapped:
       continue

    discard update_genomic_bins(record,coverage_bins,bin_size)

discard print_bed(b,bin_size,coverage_bins)
