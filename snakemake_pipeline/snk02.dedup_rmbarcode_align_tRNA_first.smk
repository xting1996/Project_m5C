
##01.samples
####--------------------------------------------------------
##5hmC & 5fC total RNA 分析
##------------Version1------------------------
##2022-05-11
##seq in Hangzhou
## add spike-in seq for quantitive
####--------------------------------------------------------

#SAMPLES=[
#"m5C_Spike_in_Ctrl",
#"m5C_Spike_in_Label",
#]


SAMPLES=[
"Ctrl_rep1",
"Ctrl_rep2",
"Label_rep1",
"Label_rep2",
]
## Main pipeline
## 1.remove adapter
## 2.remove dup
## 3.remove R2 前10nt UMI
## 4.进行align
## 5.进行tRF的分析
####--------------------------------------------------------
#CUTADAPT="/lustre2/chengqiyi_pkuhpc/zhangxt/software/cutadapt/bin/cutadapt"
TRIM_GALORE = "/path/to/trim_galore"
SEQKIT = "/path/to/seqkit"
SEQTK = "/path/to/seqtk"
hisat3n = "/path/to/hisat-3n"


##03.reference
##比对到mouse p6 + spikein 


mm10_tRNA = "/path/to/smRNA"

####--------------------------------------------------------
rule all:
    input:
        expand("../02Trim/{sample}_R2.fq.gz",sample=SAMPLES),
        expand("../03rmdup/{sample}.rmdup.R2.fq.gz",sample=SAMPLES),
        expand("../03rmdup/{sample}.R2.dups.fq.gz",sample=SAMPLES),
        expand("../04rmdup_rmUMI/{sample}.rmdup.rmUMI.R2.fq.gz",sample=SAMPLES),
        expand("../05Align_tRNA/{sample}.hisat3n.tRNA.sam",sample=SAMPLES),
        expand("../05Align_tRNA/{sample}.hisat3n.tRNA.sort.bam",sample=SAMPLES),
        expand("../05Align_tRNA/{sample}.hisat3n.tRNA.sort.bam.bai",sample=SAMPLES),
        expand("../06_derRNA_file/{sample}.detRNA.fq.gz",sample=SAMPLES),
rule rmdup:
    input:
        "../02Trim/{sample}_R2.fq.gz",
    output:
        "../03rmdup/{sample}.rmdup.R2.fq.gz",
        temp("../03rmdup/{sample}.R2.dups.fq.gz")
    shell:
        """
        {SEQKIT} rmdup -s {input} -o {output[0]} -d {output[1]}
        """
rule deRandom10N_clean:
        input:
                "../03rmdup/{sample}.rmdup.R2.fq.gz",
        output:
                "../04rmdup_rmUMI/{sample}.rmdup.rmUMI.R2.fq.gz",
        shell:
                """
                #{FASTX} -Q 33 -f 11 -i {input} -o {output}
                #fastp -i {input} -o {output} -U --umi_loc=read1 --umi_len=4
                fastp -i {input} -o {output} -f 10 -t 10 -l 20
                """
                
####----------------------Hisat2-3n------------------------------

rule align_smallRNA:
        input:
                 "../04rmdup_rmUMI/{sample}.rmdup.rmUMI.R2.fq.gz",
        output:
                temp("../05Align_tRNA/{sample}.hisat3n.tRNA.sam"),
                "../06_derRNA_file/{sample}.detRNA.fq.gz"
        log:
                "../05Align_tRNA/{sample}.hisat3n.log"
        shell:
                """
                {hisat3n} \
                -x {mm10_tRNA} \
                -U {input} \
                -S {output[0]} \
                --base-change C,T \
                --directional-mapping \
                --repeat \
                -p 20 \
                --un-gz {output[1]} 
                """
rule hisat3n_sam2bam:
        input:
                "../05Align_tRNA/{sample}.hisat3n.tRNA.sam"
        output:
                "../05Align_tRNA/{sample}.hisat3n.tRNA.sort.bam"
        shell:
                "samtools sort {input} -o {output}"

rule buildIndex:
        input:
                "../05Align_tRNA/{sample}.hisat3n.tRNA.sort.bam"
        output:
                "../05Align_tRNA/{sample}.hisat3n.tRNA.sort.bam.bai"
        shell:
                "samtools index {input}"
