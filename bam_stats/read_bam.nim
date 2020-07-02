import algorithm
import hts
import os
import math
import stats

# open a bam/cram and look for the index.
let file_name = paramStr(1)
var b:Bam
open(b, file_name)

let
 mapq=10'u
 innie=0

var i=0
var reads=10000000
var insert_size: seq[int64]
var max_size= int64(100000)

var insert_size_stats: RunningStat


for record in b:

  if record.flag.dup:
     continue
  elif record.flag.secondary:
     continue
  elif record.flag.supplementary:
     continue
  elif record.flag.unmapped:
     continue


  if record.mapping_quality > mapq and record.flag.pair and not record.flag.mate_unmapped and not (record.mate_pos < record.start or record.chrom != record.mate_chrom):
    if record.mate_pos - record.start > max_size:
      continue
    i+=1
    #echo record.chrom , " " , record.start, " " , record.stop
    insert_size.add(record.mate_pos - record.start)
    insert_size_stats.push( int(record.mate_pos - record.start) )

    if i > reads:
      break

sort(insert_size, system.cmp)


echo "Insert size median: ", insert_size[ int64(i/2)-1]
echo "Insert size 99.5 percentile: ", insert_size[ i-int64(i/995)-1]


echo "Insert size mean: ", insert_size_stats.mean()
echo "Insert size std: ", insert_size_stats.standardDeviation()
