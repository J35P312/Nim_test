import algorithm
import hts
import os


# open a bam/cram and look for the index.
let file_name = paramStr(1)
var b:Bam
open(b, file_name)

let
 mapq=10'u
 max_size=100000'u
 innie=0

var i=0
var reads=1000000
var insert_size: seq[int64]

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
    i+=1
    #echo record.chrom , " " , record.start, " " , record.stop
    insert_size.add(record.mate_pos - record.start)

    if i > reads:
      break

sort(insert_size, system.cmp)
echo insert_size[ int64(i/2)-1]
echo insert_size[ i-int64(i/99)-1]
