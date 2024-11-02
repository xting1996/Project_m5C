
SAMPLES=[
"Ctrl-rep1",
"Ctrl-rep2",
"Label-rep1",
"Label-rep2",
]
TRIM_GALORE = "/path/to/trim_galore"
####--------------------------------------------------------
rule all:
    input:
        expand("../01_rawData/{sample}_R1.fq.gz",sample=SAMPLES),
        expand("../01_rawData/{sample}_R2.fq.gz",sample=SAMPLES),
        expand("../02Trim/{sample}_R1_val_1.fq.gz",sample=SAMPLES),
        expand("../02Trim/{sample}_R2_val_2.fq.gz",sample=SAMPLES),

rule trim:
    input:
        fq1 = "../01_rawData/{sample}_R1.fq.gz",
        fq2 = "../01_rawData/{sample}_R2.fq.gz",
    output:
        "../02Trim/{sample}_R1_val_1.fq.gz",
        "../02Trim/{sample}_R2_val_2.fq.gz",
    params:
        out_path = "../02Trim/"
    shell:
        """
        {TRIM_GALORE} -q 20 --fastqc --phred33 \
        --paired {input.fq1} {input.fq2} -o {params.out_path} \
        --cores 20 --length 30  
        """
