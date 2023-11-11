####
####Xiaoting Zhang
####filter Hisat-3N alignments conversion numbers
## Yf:i:N
## Yf is the conversion number record
####-----------------Version 1--------------------------------
##   1.YF,clean conversion base numbers
##   Yf:i:N
##   Yf is the conversion number record
##   2.date @ 2023-05-24

##-----------------Version 2--------------------------------
## 1.clean other mismatches
## 2.clean total mismatches
## XM:i:<N> : The number of mismatches in the alignment. Only present if SAM record is for an aligned read.



import pysam
import argparse

##拆文件
parser = argparse.ArgumentParser(description="According full length tRNA BAM(60-150nt), to remove non-templated base in short reads aligned BAM,generate a new fastq file")


parser.add_argument("-input_raw_BAM", "--Input_the_raw_file",
                    help="Input the raw BAM file ",required=True)

parser.add_argument("-output_filter_BAM", "--output_the_filter_BAM_file",
                    help="Output the the filter BAM file ",required=True)

ARGS = parser.parse_args()



Input_file = ARGS.Input_the_raw_file
Output_file = ARGS.output_the_filter_BAM_file

inputBam = pysam.AlignmentFile(Input_file,"rb")

f_out = pysam.AlignmentFile(Output_file,"wb", template=inputBam)


for read in inputBam.fetch(until_eof=True):
    if not read.is_unmapped:
        seq_length = len(read.query_sequence)
        align_length = len(read.query_alignment_sequence)
        conversion_base = read.get_tag("Yf")
        non_conversion_base = read.get_tag("Zf")
        total_mismatch = read.get_tag("XM")
#        able_conversion = conversion_base + non_conversion_base + 1
#        conversion_ratio_in_reads = conversion_base*1.0 / align_length

#        if ((align_length >= 20) and (conversion_base <= 4) and (conversion_ratio_in_reads <= 0.1)):
        #if ((align_length >= 20) and (conversion_base <= 3) and (total_mismatch <= 3) ):
        if ((align_length >= 20) and (conversion_base <= 2) and (total_mismatch <= 2) ):
            f_out.write(read)

inputBam.close()
f_out.close()
