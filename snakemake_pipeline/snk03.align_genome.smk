
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
hisat3n = "/path/to/hisat-3n"
##03.reference

genome_ref = "/path/to/genome"
####--------------------------------------------------------
rule all:
    input:
        expand("../06_derRNA_file/{sample}.detRNA.fq.gz",sample=SAMPLES),
        expand("../07Align_genome/{sample}.hisat3n.genome.sam",sample=SAMPLES),
        expand("../07Align_genome/{sample}.hisat3n.genome.sort.bam",sample=SAMPLES),
        expand("../07Align_genome/{sample}.hisat3n.genome.sort.bam.bai",sample=SAMPLES),





                
####----------------------Hisat2-3n------------------------------

rule align_HISAT3N_genome:
        input:
                "../06_derRNA_file/{sample}.detRNA.fq.gz"
        output:
                temp("../07Align_genome/{sample}.hisat3n.genome.sam"),
        shell:
                """
                {hisat3n} \
                -x {mouse_ref} \
                -U {input} \
                -S {output[0]} \
                --base-change C,T \
                --directional-mapping \
                --repeat \
                -p 20 \
                -k 1 \
                """
rule hisat3n_sam2bam:
        input:
                "../07Align_genome/{sample}.hisat3n.genome.sam"
        output:
                "../07Align_genome/{sample}.hisat3n.genome.sort.bam"
        shell:
                "samtools sort {input} -o {output}"

rule buildIndex:
        input:
                "../07Align_genome/{sample}.hisat3n.genome.sort.bam"
        output:
                "../07Align_genome/{sample}.hisat3n.genome.sort.bam.bai"
        shell:
                "samtools index {input}"
