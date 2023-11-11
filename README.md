# m5C-TAC-tools V1.0

### What is m5C-TAC-seq?

m5C-TAC-seq is a bisulfite-free approach that combines TET-assisted m5C-to-f5C oxidation with selective chemical labeling, enabling pre-enrichment and C-to-T transitions at m5C sites, this approach is mild on RNA and does not affect unmodified Cs, thus allowing direct m5C detection in low-abundance and low-sequence-complexity RNAs. m5C-TAC-seq is a sensitive, accurate, and robust method for transcriptome-wide m5C detection.

### What is m5C-TAC-tool ?
The m5C-TAC-seq tool can accurately identify m5C sites from m5C-TAC-seq data generated from experiments.

<img width="1000" alt="image" src="https://github.com/xting1996/Project_m5C/assets/34152806/1a754c02-f3a6-4a36-815b-80a7c24923a9">


### Update log
1. 2023-11-11，first commit，by Xiaoting Zhang

##  Background

### 1.m5C-TAC-seq and m5C-TAC-tools

m5C-TAC-seq can direct detect m5C sites through TET2-mediated m5C-to-f5C oxidation and AI-labeling， the m5C sites are read as "T" in NGS data， thus m5C-TAC-tools parse the  C-to-T signal between treated-sample and control-sample to call the accurate m5C sites.



### 2. Pre-processing  the environment and software

1. TrimGalore (v0.6.5)
2. fastq-multx (v1.4.3)
3. seqkit (v0.10.0)
4. fastp (v0.23.2)
5. HISTA3N (v2.2.1-3n-0.0.3)
6. pysamstats (v0.14)
7. snakemake (v7.30.1)
8. Python (v 3.8.13) and pysam (v 0.16.0.1).
### 3. Pre-processing the raw sequencing data

m5C-TAC-seq is based on a multiplexing library construction strategy, therefore reducing the input of starting material and minimizing the technical variations among different samples.

If you didn't use the multiplexing library construction strategy, please ignore it, which is the following step2.

##### step1. trim adapter
```
TRIM_GALORE = /path/to/trim_glore
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
        
```
`
##### step2. Assign merged reads to respective individual samples according to sample-barcoded sequence.
```
fastq-multx -B barcode.txt -m 0 \
mergesample_R1_val_1.fq.gz mergesample_R2_val_2.fq.gz \
-o %_R1.fq.gz -o %_R2.fq.gz
```

the `barcode.txt` file records the sample-barcoded information with following format:
```
sample1    CATCTG
sample2    AGGTAC
sample3    TGTGCT
.....      ......
samplen    NNNNNN
```

##### step3.remove duplication reads

```
SEQKIT = /path/to/seqkit
rule rmdup:
    input:
        "../02Trim/{sample}_R2_val_2.fq.gz",
    output:
        "../03rmdup/{sample}.rmdup.R2.fq.gz",
    shell:
        """
        {SEQKIT} rmdup -s {input} -o {output[0]} 
        """
        
```

##### step4 .remove UMI
```
FASTP = /path/to/fastp
rule deRandom10N_clean:
        input:
                "../03rmdup/{sample}.rmdup.R2.fq.gz",
        output:
                "../04rmdup_rmUMI/{sample}.rmdup.rmUMI.R2.fq.gz",
        shell:
                """
                {FASTP} -i {input} -o {output} -f 10 -t 10
                """
```

### 4. Pre-processing  the references 

###### 1. small RNA reference
To identify the transcriptome-wide m5C sites and tRNA m5C sites, we firstly align the reads to small RNA references，including tRNA，rRNA, snRNA, snoRNA, pre-miRNA，mature-miRNA, etc.

```
## genomic tRNA,download from GtRNAdb,hg19
## pre-tRNA,need add "CCA" to the end of each tRNA seq
wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-mature-tRNAs.fa

## mito-tRNA,download from mitotRNAdb
```

rRNA/snRNA/snoRNA are download from NCBI by searching the key words. and pre-miRNA/mature-miRNA are download from miRBase21.

Different small RNA references are merged into `smallRNA.fa`
###### 2.  genome reference
+ genome file 
```
## human,hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

## mouse,mm10

wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
```

+ gtf file

```
## human, hg19 gtf
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

## mouse, mm10 gtf

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz

```

## m5C-TAC-tools for m5C-TAC-seq analysis

#### step 01. Download reference genome and build index
```
## for small RNA
hisat-3n-build \
    --base-change C,T \
    -p 20 \
    smallRNA.fa smallRNA


## for genome
hisat-3n-build \
    --base-change C,T \
    -p 20 \
    --repeat-index 50-300  genome.fa genome_repeat

``` 


#### step02. Map the m5C-TAC-seq data

##### 1. Removing the reads which aligned to small RNA
```
hista3n = /path/to/hisat-3n

rule align_smallRMA:
        input:
                 "../04rmdup_rmUMI/{sample}.rmdup.rmUMI.R2.fq.gz",
        output:
                temp("../05Align_tRNA/{sample}.hisat3n.tRNA.sam"),
                "../06_derRNA_file/{sample}.detRNA.fq.gz"
        log:
                "../05Align_tRNA/{sample}.hisat3n.log"
        shell:
                """
                {hista3n} \
                -x {mm10_tRNA} \
                -U {input} \
                -S {output[0]} \
                --base-change C,T \
                --directional-mapping \
                --repeat \
                -p 20 \
                --un-gz {output[1]}
                """
```

##### 2. Align to genome

```
hista3n = /path/to/hisat-3n

rule align_HISAT3N_genome:
        input:
                "../06_derRNA_file/{sample}.detRNA.fq.gz"
        output:
                temp("../07Align_genome/{sample}.hisat3n.genome.sam"),
        shell:
                """
                {hista3n} \
                -x {mouse_ref} \
                -U {input} \
                -S {output[0]} \
                --base-change C,T \
                --directional-mapping \
                --repeat \
                -p 20 \
                -k 1 \
                """
                
```

##### step03. bam filter
```
rule bam_filter:
        input:
                "../07Align_genome/{sample}.sort.bam"
        output:
                "../07Align_genome/{sample}.filter.mis2.sam",
        shell:
                """
                python 01_clean_Hisat3n_m5C_align.py \
                -input_raw_BAM {input} \
                -output_filter_BAM {output}
                """
                

```

#### step03. Parsing the bam file

##### step01. parse the bam file by pysamstats
```
rule bam2pmat:
    input:
        bai = "../07Align_genome/{sample}.filter.mis2.{strand}.bam.bai",
        bam = "../07Align_genome/{sample}.filter.mis2.{strand}.bam",
    output:
        "../07Parse/{sample}.{RST}.{strand}.pmat"
    shell:
        """
         pysamstats -D 100000 -S "nofilter" -t variation {input.bam} --fasta {HG19_FA} --chromosome {wildcards.RST} > {output}
        """
```

##### step02. tidy the output .pmat file to get the C position
```
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
```

#### step03. Calling m5C sites

##### step01. merge the control and treated sample
```
Rscript script03_pmat_merge_ctrl_and_label.r
```

##### step02. Set cutoff for m5C calling

This step is based on a custom cutoff.

### Maintainers and Contributing

m5C-TAC-tools is developed and maintained by Xiaoting Zhang (xiaotingzhang@pku.edu.cn).
