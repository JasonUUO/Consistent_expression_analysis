SRA,FRR = glob_wildcards("rawReads/{sra}_{frr}.fastq.gz")

rule all:
    input:
        expand("rawQC/{sra}_{frr}_fastqc.{extension}", sra=SRA, frr=FRR,extension=["zip","html"]),
        expand("multiqc_report.html"),
        expand("trimmedreads{sra}_fastq.html", sra=SRA),
        expand("genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"),
        expand("genome/Mus_musculus.GRCm39.106.gtf"),
        expand("starAlign/{sra}Aligned.sortedByCoord.out.bam", sra=SRA),
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
        fazip = "genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
    shell:
        "cd genome"
        " && wget ftp://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"

rule bwa_index:
    "Generating index with BWA"
    input:
        fazip = rules.get_genome_fa.output.fazip
    output:
        #fa = "genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa",
        amb = "genome/ref.amb", 
        ann = "genome/ref.ann", 
        bwt = "genome/ref.bwt", 
        pac = "genome/ref.pac", 
        sa = "genome/ref.sa"
    shell:
        "gunzip -k {input.fazip} "
        "&& bwa index -a bwtsw {input.fazip} -p genome/ref "
        "&& rm genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa" 
    

rule mem_bwa:
    "Aligning with BWA MEM"
    input:
        index = rules.bwa_index.output,
        fastq1 = rules.fastp.output.read1,
        fastq2 = rules.fastp.output.read2
    output:
        bam = "aligned/{sra}.bam"
    threads: 
        32
    shell:
        "bwa mem -t {threads} genome/ref {input.fastq1} {input.fastq2} | samtools sort -@ {threads} -o {output.bam}"

        
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
        
