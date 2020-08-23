import hts
import sequtils
import tables
import math
import system
import strutils

proc make_coverage_structure(b: Bam, bin_size: int) :(Table[string,seq[float]],Table[string,seq[float]],Table[string,seq[float]],Table[string,seq[float]] )=
  var coverage_bins = initTable[string, seq[float]]()
  var spanning_read_bins = initTable[string, seq[float]]()
  var discordant_pair_bins = initTable[string, seq[float]]()
  var quality_bins = initTable[string, seq[float]]()

  var header= b.hdr
  for line in targets(header):
    coverage_bins[line.name]= newSeqWith( int(ceil( float(line.length)/float(bin_size))) , float(0))
    spanning_read_bins[line.name]= newSeqWith( int(ceil( float(line.length)/float(bin_size))) , float(0))
    discordant_pair_bins[line.name]= newSeqWith( int(ceil( float(line.length)/float(bin_size))) , float(0))
    quality_bins[line.name]= newSeqWith( int(ceil( float(line.length)/float(bin_size))) , float(0))

  return(coverage_bins,spanning_read_bins,discordant_pair_bins,quality_bins)

proc update_genomic_bins(alignment: Record,
  genome_bins:var Table[string,seq[float]],
  spanning_read_bins:var Table[string,seq[float]],
  discordant_pair_bins:var Table[string,seq[float]],
  quality_bins:var Table[string,seq[float]],
  bin_size:int): int=

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

proc print_bed(b:Bam,
  bin_size:int,
  genome_bins:var Table[string,seq[float]],
  spanning_read_bins:var Table[string,seq[float]],
  discordant_pair_bins:var Table[string,seq[float]],
  quality_bins:var Table[string,seq[float]]): string =

  var output_bed="#chromosome\tstart\tend\tcoverage\n"
  for chrom in targets(b.hdr):
    for i in 0 .. genome_bins[chrom.name].len-1:
      output_bed.add( chrom.name & "\t" & $(1+i*bin_size) & "\t" & $((i+1)*bin_size) & "\t" & $genome_bins[chrom.name][i] & "\n")

  return(output_bed.strip(leading = false))

proc print_wig(b:Bam,
  bin_size:int,
  genome_bins:var Table[string,seq[float]],
  spanning_read_bins:var Table[string,seq[float]],
  discordant_pair_bins:var Table[string,seq[float]],
  quality_bins:var Table[string,seq[float]]): string =

  var output_wig="track type=wiggle_0 name=\"Coverage\" description=\"Per bin average coverage\"\n"
  for chrom in targets(b.hdr):
    add(output_wig,"fixedStep chrom=" & chrom.name & " start=1 step=" & $bin_size & "\n")
    for i in 0 .. genome_bins[chrom.name].len-1:
      output_wig.add( $genome_bins[chrom.name][i] & "\n")

  return(output_wig.strip(leading = false))

export print_bed 
export print_wig 
export make_coverage_structure 
export update_genomic_bins
