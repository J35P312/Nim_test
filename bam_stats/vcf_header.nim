import hts/bam
import os
import strutils

proc find_sample(bam_header: Header): string =

 var txt = newString(bam_header.hdr.l_text)
 copyMem(txt[0].addr, bam_header.hdr.text, txt.len)

 for line in txt.split("\n"):  
  if line.startsWith("@RG") and "\tSM:" in line:
   return( line.split("\tSM:")[1].split("\t")[0].strip() )

 return("Sample")

proc print_header(libraryData:string, version: string,bam_header:Header): string =

 var headerString="##fileformat=VCFv4.1\n"
 headerString.add("##source=TIDDIT-" & version &  "\n")

 #define the alowed events
 headerString.add("##ALT=<ID=DEL,Description=\"Deletion\">\n")
 headerString.add("##ALT=<ID=DUP,Description=\"Duplication\">\n")
 headerString.add("##ALT=<ID=TDUP,Description=\"Tandem duplication\">\n")
 headerString.add("##ALT=<ID=INV,Description=\"Inversion\">\n")
 headerString.add("##ALT=<ID=INS,Description=\"Insertion\">\n")
 headerString.add("##ALT=<ID=BND,Description=\"Break end\">\n")

 for contig in targets(bam_header):
   headerString.add("##contig=<ID=" & contig.name & ",length=" & $contig.length & ">\n")

 #Define the info fields
 headerString.add("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
 headerString.add("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of an intra-chromosomal variant\">\n")
 headerString.add("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
 headerString.add("##INFO=<ID=LFA,Number=1,Type=Integer,Description=\"Links from window A\">\n")
 headerString.add("##INFO=<ID=LFB,Number=1,Type=Integer,Description=\"Links from window B\">\n")
 headerString.add("##INFO=<ID=LTE,Number=1,Type=Integer,Description=\"Links to event\">\n")
 headerString.add("##INFO=<ID=COVA,Number=1,Type=Float,Description=\"Coverage on window A\">\n")
 headerString.add("##INFO=<ID=COVM,Number=1,Type=Float,Description=\"The coverage between A and B\">\n")
 headerString.add("##INFO=<ID=COVB,Number=1,Type=Float,Description=\"Coverage on window B\">\n")
 headerString.add("##INFO=<ID=OR,Number=4,Type=Integer,Description=\"Orientation of the pairs (FF,RR,RF,FR)\">\n")
 headerString.add("##INFO=<ID=ORSR,Number=2,Type=Integer,Description=\"Orientation of the split reads (inverted,normal)\">\n")
 headerString.add("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
 headerString.add("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n")
 headerString.add("##INFO=<ID=QUALA,Number=1,Type=Float,Description=\"The average mapping quality of the reads in window A\">\n")
 headerString.add("##INFO=<ID=QUALB,Number=1,Type=Float,Description=\"The average mapping quality of the reads in window B\">\n")

 #set filters
 headerString.add("##FILTER=<ID=BelowExpectedLinks,Description=\"The number of links or reads between A and B is too small\">\n")
 headerString.add("##FILTER=<ID=FewLinks,Description=\"Unexpectedly low fraction of discordant reads betwen A and B\">\n")
 headerString.add("##FILTER=<ID=UnexpectedCoverage,Description=\"The coverage of the window on chromosome B or A is higher than 4*average coverage\">\n")
 headerString.add("##FILTER=<ID=Smear,Description=\"Window A and Window B overlap\">\n")
 headerString.add("##FILTER=<ID=RegionalQ,Description=\"The mapping quality of the region is lower than the user set limit\">\n")
 headerString.add("##FILTER=<ID=MinSize,Description=\"The variant is smaller than the user set limit\">\n")
 headerString.add("##FILTER=<ID=Ploidy,Description=\"Intrachromosomal variant on a chromosome having 0 ploidy\">\n")
 headerString.add("##FILTER=<ID=SplitsVSDiscs,Description=\"large variant supported mainly by split reads (and not discorant pairs) \">\n")
 headerString.add("##FILTER=<ID=Density,Description=\"The discordant reads cluster too tightly\">\n")

 #set format 
 headerString=headerString & "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
 headerString=headerString & "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n"
 headerString=headerString & "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of paired-ends that support the event\">\n"
 headerString=headerString & "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"Number of split reads that support the event\">\n"
 headerString=headerString & "##FORMAT=<ID=DR,Number=2,Type=Integer,Description=\"Number of paired-ends that supporting the reference allele (breakpoint A, and B)\">\n"
 headerString=headerString & "##FORMAT=<ID=RR,Number=2,Type=Integer,Description=\"Number of reads supporting the reference allele (breakpoint A, and B)\">\n"

 #this string contains info such as insert size and mean coverage
 headerString.add(libraryData & "\n")

 headerString.add("TIDDITcmd=\"" & getAppFilename() & " " & commandLineParams().join() & "\"\n")
 headerString.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" & find_sample(bam_header) & "\n")
 return(headerString)

export print_header
