SRA,FRR = glob_wildcards("rawReads/{sra}_{frr}.fastq.gz")

rule all:
    input:
        expand("rawQC/{sra}_{frr}_fastqc.{extension}", sra=SRA, frr=FRR,extension=["zip","html"]),
        expand("multiqc_report.html"),
        expand("trimmedreads{sra}_fastq.html", sra=SRA),
        "genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa",
        "genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa",
        expand("aligned/{sra}.bam", sra=SRA),
        expand("logs/{sra}_sum.txt", sra=SRA),
        expand("logs/{sra}_met.txt", sra=SRA),
        ["index."  + str(i) + ".ht2" for i in range(1,9)]
        expand("rawcounts/rawcounts.tsv",)
        expand("AML_gene_lists.csv",)
        
        
rule rawFastqc:
    input:
        rawread="rawReads/{sra}_{frr}.fastq.gz",
    output:
        zip="rawQC/{sra}_{frr}_fastqc.zip",
        html="rawQC/{sra}_{frr}_fastqc.html",
    threads:
        1
    params:
        path="rawQC/",
    shell:
        """
        fastqc {input.rawread} --threads {threads} -o {params.path}
        """
        
rule multiqc:
    input:
        rawqc=expand("rawQC/{sra}_{frr}_fastqc.zip",sra=SRA ,frr=FRR),
    output:
       "multiqc_report.html"
    shell:
        """
        multiqc {input.rawQC}
        """
        
rule fastp:
     input:
         read1="rawReads/{sra}_1.fastq.gz",
         read2="rawReads/{sra}_2.fastq.gz",
     output:
         read1="trimmedreads/{sra}_1P.fastq.gz",
         read2="trimmedreads/{sra}_2P.fastq.gz",
         report_html= "trimmedreads{sra}_fastq.html",
     threads: 
        4
     shell:
         """
         fastp --thread {threads} -i {input.read1} -I {input.read2} -o {output.read1} -O {output.read2} -h {output.report_html}
         """
      
      
rule get_genome_fa:
    "Downloading Genome sequence, Mus Musculus primary assembly (GRCm39)"
    output:
        fazip = "genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz",
        fa = "genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
    shell:
        "cd genome"
        " && wget ftp://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
        " && gunzip -k {output.fazip} "

rule index:
    input:
        fa = rules.get_genome_fa.output.fa
    output:
        dir = ["index."  + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    threads:
        16
    shell:
        " hisat2-build -p {threads} {input.fa} index --quiet"
        " && rm genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa" 
    

rule hisat_align:
    input:
        fastq1 = rules.fastp.output.read1,
        fastq2 = rules.fastp.output.read2,
        index = rules.index.output.dir
    output:
        bams  = "aligned/{sra}.bam",
        sum   = "logs/{sra}_sum.txt",
        met   = "logs/{sra}_met.txt"
    message:
        "mapping reads to genome to bam files."
    threads: 
        16
    shell:
        "hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x index \
        -1 {input.fastq1} -2 {input.fastq2} | samtools view -Sb -F 4 -o {output.bams}"

        
rule bamtools_filter_json:
    input:
        "{sample}.bam"
    output:
        "filtered/{sample}.bam"
    params:
        json="filtering-rules.json",
        region="" # optional parameter for defining a specific region, e.g. "chr1:500..chr3:750"
    log:
        "logs/bamtools/filtered/{sample}.log"
    wrapper:
        "v1.5.0/bio/bamtools/filter_json"

        
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

region = snakemake.params.get("region")
region_param = ""

if region and region is not None:
    region_param = ' -region "' + region + '"'
        
rule htseq_count:
    input: 
        alignment = rules.star_align.output.bam 
        gtf_1 = rules.get_genome_gtf.output.gtf 
    output:
        raw_count = "rawcounts/rawcounts.csv"
    log:
        out = "rawcounts/rawcounts.out"
        err = "rawcounts/rawcounts.err"
    shell:
        "htseq-count {input.alignment} {input.gtf_1} -i gene_id --add-chromosome-info -m union"

rule normalise_and_DE_analysis:
    input:
        rules.htseq_count.output.raw_count
    output:
        "AML_gene_lists.csv"
    script:
        "Analysis.Rmd"
        
