import hts/bam
import sequtils
import tables

proc generate_signal_structure(bam_header: Header): (Table[string,seq[seq[uint64]]], Table[string,seq[string]], Table[string,seq[seq[uint64]]], Table[string,seq[string]])=

  var discordant_pairs:Table[string,seq[seq[uint64]]]
  var discordant_pairs_read_name:Table[string,seq[string]]
  var split_reads:Table[string,seq[seq[uint64]]]
  var split_reads_read_name:Table[string,seq[string]]
  
  for chromA in targets(bam_header):
    for chromB in targets(bam_header):
      if chromB.name < chromA.name:
        continue

      discordant_pairs[chromA.name & "_vs_" & chromB.name ]= @[]
      discordant_pairs_read_name[chromA.name & "_vs_" & chromB.name ]= @[]
      split_reads[chromA.name & "_vs_" & chromB.name ]= @[]
      split_reads_read_name[chromA.name & "_vs_" & chromB.name ]= @[]
 
  return(discordant_pairs, discordant_pairs_read_name, split_reads, split_reads_read_name)

proc sv_signals(alignment: Record,max_insert_size:int64,mapq:uint8,
  discordant_pairs:var Table[string,seq[seq[uint64]]],
  discordant_pairs_read_name:var Table[string,seq[string]],
  split_reads:var Table[string,seq[seq[uint64]]],
  split_reads_read_name:var Table[string,seq[string]]): string =

  var signal=""

  if alignment.mapping_quality < mapq:
    return("")

  if alignment.flag.mate_unmapped or (not alignment.flag.pair) or alignment.flag.supplementary or alignment.flag.secondary:
    return("")

  var orientation:string = "FF"
  if (alignment.flag.reverse and alignment.flag.mate_reverse):
    orientation="RR"
  elif (alignment.flag.reverse and not alignment.flag.mate_reverse):
    orientation="RF"
  elif (not alignment.flag.reverse and alignment.flag.mate_reverse):
    orientation="FR"

  if (alignment.chrom < alignment.mate_chrom):
    signal.add( alignment.chrom & "\t" & $alignment.start & "\t" & $alignment.stop & "\t" & alignment.mate_chrom & "\t" & $alignment.mate_pos & "\t" & orientation & "\n")
  elif (alignment.chrom == alignment.mate_chrom and int64(alignment.mate_pos-alignment.start) > max_insert_size):
    signal.add( alignment.chrom & "\t" & $alignment.start & "\t" & $alignment.stop & "\t" & alignment.mate_chrom & "\t" & $alignment.mate_pos & "\t" & orientation & "\n")
  return(signal)

export generate_signal_structure, sv_signals
