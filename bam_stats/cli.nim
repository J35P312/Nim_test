import read_bam
import vcf_header
import hts
import os

var output_prefix = "output"
let version = "0.0.0"

let (median_ins,ins_treshold,library_type,avg_read_length)= bam_stats(paramStr(1),50000000)

var file_name:string=paramStr(1)
var b:Bam
open(b, file_name)

let library_stats="##LibraryStats=TIDDIT-" & version & " Coverage=0.0" & " AverageReadLength=" & $avg_read_length & " MedianInsertSize=" & $median_ins & " 99.5PercentileInsertSize=" & $ins_treshold & " LibraryType=" & library_type
#let library_stats="##LibraryStats=TIDDIT-" & version & " Coverage=0.0"

let output_vcf_header = print_header(library_stats,version,b.hdr)
writeFile( output_prefix & ".vcf", output_vcf_header)
