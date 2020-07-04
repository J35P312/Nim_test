import algorithm
import hts
import os
import math
import stats


proc bam_stats(file_name: string): int =
  var b:Bam
  open(b, file_name)

  let
    mapq=10'u


  var analysed_pairs=0

  var i=0
  var max_size = 100000
  var reads=50000000
  #var reads=500000

  var insert_size: seq[int64]
  var innie = uint64(0)
  var outtie = uint64(0)
  var pairs = uint64(0)
  var sequence_str: string

  var insert_size_stats: RunningStat
  var read_length_stats: RunningStat


  for record in b:

    if record.flag.dup:
       continue
    elif record.flag.secondary:
       continue
    elif record.flag.supplementary:
       continue
    elif record.flag.unmapped:
       continue

    if record.flag.pair and record.flag.read1:
      pairs+=1

    if record.mapping_quality > mapq and record.flag.pair and not record.flag.mate_unmapped and not (record.mate_pos < record.start or record.chrom != record.mate_chrom):

      if record.mate_pos - record.start > max_size:
        continue

      insert_size.add(record.mate_pos - record.start)
      insert_size_stats.push( int(record.mate_pos - record.start) )
      analysed_pairs+=1

      if record.flag.reverse and not record.flag.mate_reverse:
        outtie+=1
      elif record.flag.mate_reverse and not record.flag.reverse:
        innie+=1

    read_length_stats.push( record.sequence(sequence_str).len )
    i+=1
    if i > reads:
      break

  sort(insert_size, system.cmp)

  echo "sampled read-pairs:", pairs
  echo "Insert size median: ", insert_size[ int64(analysed_pairs/2)-1]
  echo "Insert size 99.9 percentile: ", insert_size[ analysed_pairs-int64(analysed_pairs/999)-1]
  echo "Insert size mean: ", insert_size_stats.mean()
  echo "Insert size std: ", insert_size_stats.standardDeviation()
  echo ""
  echo "Pair orientation --> <--: " , innie
  echo "Pair orientation <-- -->: " , outtie
  if innie > uint64(0) or outtie > uint64(0):
    if innie > outtie:
      echo "major orientation: --> <--"
    else:
      echo "major orientation: <-- -->"
  else:
    echo "Single-end"

  echo ""
  echo "sampled reads:", i
  echo "average read length:", read_length_stats.mean()
  echo "read length standard deviation:", read_length_stats.standardDeviation()
  echo "========================="

export bam_stats

var slask: int
slask= bam_stats(paramStr(1))
