# _*_ coding: UTF-8 _*_

########################################################################
# MENG Howard
# 2021-07-08
# ABE data 
######################################################################## 
# run on abyss
# /home/menghaowei/menghw_HD/ABE_project/test.01.hisat3n_ABE_data


# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>


# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
HG38_FA = "/lustre2/chengqiyi_pkuhpc/zhangxt/genome/mouse/mm10_m5C_spikein/mm10.p6.spikein.fa"
# rm INDEL info
#ABE_SNP_MOCK_INFO = "/lustre2/chengqiyi_pkuhpc/zhangxt/ABE/20220510_direct_seq/script/293T_snp/293T-Mock-Input_A_or_T.mut.UnSort.pmat.AG_or_TC.out.bed"

# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES=[
"Ctrl_rep1",
"Ctrl_rep2",
"Label_rep1",
"Label_rep2",

]
CHROM = [
"chr1",
"chr2",
"chr3",
"chr4",
"chr5",
"chr6",
"chr7",
"chr8",
"chr9",
"chr10",
"chr11",
"chr12",
"chr13",
"chr14",
"chr15",
"chr16",
"chr17",
"chr18",
"chr19",
"chrX",
"chrY",
"chrM",]

STRAND=['rev','fwd']

rule all:
    input:
        expand("../07Align_genome/{sample}.hisat3n.genome.sort.filter.mis2.bam",sample=SAMPLES,strand = STRAND),
        expand("../07Align_genome/{sample}.filter.mis2.rev.bam",sample=SAMPLES),
        expand("../07Align_genome/{sample}.filter.mis2.fwd.bam",sample=SAMPLES),
        expand("../07Align_genome/{sample}.filter.mis2.{strand}.bam",sample=SAMPLES,strand = STRAND),
        expand("../07Align_genome/{sample}.filter.mis2.{strand}.bam.bai",sample=SAMPLES,strand = STRAND),
        expand("../08Parse_C2T/{sample}.{RST}.{strand}.mpmat",sample=SAMPLES,RST=CHROM,strand = STRAND),
        expand("../08Parse_C2T/{sample}.{RST}.{strand}.tidy.mpmat",sample=SAMPLES,RST=CHROM,strand = STRAND),

rule rev_separate_strand:
    input:
        "../07Align_genome/{sample}.hisat3n.genome.sort.filter.mis2.bam"
    output:
        "../07Align_genome/{sample}.filter.mis2.rev.bam"
    shell:
        """
        samtools view -b -f 16 {input} -o {output}
        """
rule fwd_separate_strand:
    input:
        "../07Align_genome/{sample}.hisat3n.genome.sort.filter.mis2.bam"
    output:
        "../07Align_genome/{sample}.filter.mis2.fwd.bam"
    shell:
        """
        samtools view -b -F 16 {input} -o {output}
        """
rule buildindex:
    input:
        "../07Align_genome/{sample}.filter.mis2.{strand}.bam"
    output:
        "../07Align_genome/{sample}.filter.mis2.{strand}.bam.bai"
    shell:
        "samtools index {input}"
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# mpileup 
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule bam_mpileup:
    input:
        bai = "../07Align_genome/{sample}.filter.mis2.{strand}.bam.bai",
        bam = "../07Align_genome/{sample}.filter.mis2.{strand}.bam",
    output:
        "../08Parse_C2T/{sample}.{RST}.{strand}.pmat"
    shell:
        """
         pysamstats -D 100000 -S "nofilter" -t variation {input.bam} --fasta {HG38_FA} --chromosome {wildcards.RST} > {output}
        """

rule mp_filter:
    input:
        "../08Parse_C2T/{sample}.{RST}.{strand}.pmat"
    output:
        "../08Parse_C2T/{sample}.{RST}.{strand}.tidy.pmat"
    params:
        std = "{strand}"
    shell:
        """
         python script02.tidy_pmat_v5.py -input_mpmat_file {input} -output_mpmat_file {output} -strand {params.std}
        """
