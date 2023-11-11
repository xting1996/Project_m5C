######---------------------------------------------------------------
##v3,2023-05-02
##change info
##add pavlue info


######---------------------------------------------------------------

##Binomimal Test
from scipy.stats import binom_test
"""
 coverage|mutation number
    +---+---+    a is the # treat group - mutation
    | n | k |    b is the # treat group - Not mutation
    +---+---+   
"""
####set p-value
#pi = 1/1000
#pi = 2/1000
pi = 1/1000
##这个pi要去计算
##
def BinomimalTest(k,n,pi):
    try:
        BTpvalue = binom_test(k,n,pi, alternative='greater')
        return(BTpvalue)
    except:
        print("BTpvalue-ERROR")



##parse parameters
import argparse

parser = argparse.ArgumentParser(description="ABE BAM to pmat")


parser.add_argument("-input_mpmat_file", "--Input_the_raw_mpmat_file",
                    help="Input the raw mpmat file",required=True)


parser.add_argument("-output_mpmat_file", "--output_the_filter_file",
                    help="Output the filter mpmat file", type=str)


parser.add_argument("-strand", "--input_strand_info",
                    help="input the strand info, fwd or rev ", type=str)


ARGS = parser.parse_args()

input_mpmat = ARGS.Input_the_raw_mpmat_file


input_strand = ARGS.input_strand_info

output_mpmat = ARGS.output_the_filter_file

output_mpmat = ARGS.output_the_filter_file 

output_mpmat2 = output_mpmat + ".simple"
# %%
fout = open(output_mpmat,"w")
fout2 = open(output_mpmat2,"w")
# out_header = ["chrom",	"pos",	"ref",	"reads_all",	"matches",	"mismatches",	"deletions",	"A",	"C",	"T",	"G",	"reads_all_fwd",	"matches_fwd",	"mismatches_fwd",	"deletions_fwd",	"A_fwd",	"C_fwd",	"T_fwd",	"G_fwd",	"reads_all_rev",	"matches_rev",	"mismatches_rev",	"deletions_rev",	"A_rev",	"C_rev",	"T_rev",	"G_rev","C2T_del_ratio"]
# fout.write("\t".join(out_header) + "\n")

i = 0
with open(input_mpmat,"r") as f:
    # next(f)
    for lines in f:
        line = lines.strip().split("\t")
        """
        chrom	pos	ref	reads_all	reads_pp	matches	matches_pp	mismatches	mismatches_pp	deletions	deletions_pp	insertions	insertions_pp	A	A_pp	C	C_pp	T	T_pp	G	G_pp	N	N_pp
        """
        i += 1
        if i < 2:
            outline = [line[0], line[1], line[2] ,line[3] , line[5] , line[7] , line[9] , line[13] , line[15] , line[17] , line[19] ,"mismatch_ratio","C2T_ratio","pvalue"]
            fout.write("\t".join(outline) + "\n")
            
            outline2 = [line[0], line[1], line[2] ,line[3] ,"C2T_reads","C2T_ratio","Pvalue","mismatch_reads","mismatch_ratio"]
            
            fout2.write("\t".join(outline2) + "\n")
            
            
        else:
            ref_base = line[2]
            ##reads all = matches + mismatch 
            line[3] = str(int(line[5]) + int(line[7])) 

            if input_strand == "fwd":
                ##query C2T + deltion info
                if ref_base == "C":
                    all_reads = int(line[5]) + int(line[7]) 
                    if all_reads > 0:

                        C2T_del_reads = int(line[17]) 

                        C2T_del_ratio = (int(line[17]) + 0.0)/ (all_reads +0.0) * 100

                        C2T_pvalue = BinomimalTest(C2T_del_reads,all_reads,pi)
                        mismatch_ratio = (int(line[7]) + 0.0) /(int(line[3]) + 0.0) * 100
                        ##change all reads info
                        #line[3] = str((int(line[5]) + int(line[7])))
                        
                        outline = [line[0], line[1], line[2] ,line[3] , line[5] , line[7] , line[9] , line[13] , line[15] , line[17] , line[19]]

                        newline = outline + [str(C2T_del_ratio),str(C2T_pvalue)]
                        
                        newline_simple = outline[:4] + [str(C2T_del_reads),str(C2T_del_ratio),str(C2T_pvalue),str(line[7]),str(mismatch_ratio)]

                        fout.write("\t".join(newline) + "\n")
                        fout2.write("\t".join(newline_simple) + "\n")

            elif input_strand == "rev":
                ##query G2A info
                    
                if ref_base == "G":
                    
                    all_reads = int(line[5]) + int(line[7]) 
                    
                    if all_reads > 0:
                        
                        
                        C2T_del_reads = int(line[13])

                        C2T_del_ratio = (int(line[13]) + 0.0)/ ( all_reads +0.0) * 100

                        C2T_pvalue = BinomimalTest(C2T_del_reads,all_reads,pi)
                        mismatch_ratio = (int(line[7]) + 0.0) /(int(line[3]) + 0.0) * 100
                          
                        #line[3] = str((int(line[5]) + int(line[7])))

                        outline = [line[0], line[1], line[2] ,line[3] , line[5] , line[7] , line[9] , line[13] , line[15] , line[17] , line[19]]

                        newline = outline + [str(mismatch_ratio),str(C2T_del_ratio),str(C2T_pvalue)]
                        
                        newline_simple = outline[:4] + [str(C2T_del_reads),str(C2T_del_ratio),str(C2T_pvalue),str(line[7]),str(mismatch_ratio)]
                        

                        fout.write("\t".join(newline) + "\n")
                        fout2.write("\t".join(newline_simple) + "\n")


fout.close()
fout2.close()
