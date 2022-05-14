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
    "Download Genome sequence, Mus Musculus primary assembly (GRCm39)"
    output:
        fa = "genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
    shell:
        "cd genome"
        " && wget ftp://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz "
        " && gunzip -k Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"


rule get_genome_gtf:
    "Download gtf annotations corresponding to our genome"
    output:
        gtf = "genome/Mus_musculus.GRCm39.106.gtf"
    shell:
        "cd genome && " 
        " wget ftp://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz "
        " && gunzip -k Mus_musculus.GRCm39.106.gtf.gz"


rule star_index:
    input:
        gtf = rules.get_genome_gtf.output.gtf ,
        fa = rules.get_genome_fa.output.fa       
    output:
        dir = "starIndex",
    threads:
        32
    shell:
        " STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.dir} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} "
        " && rm genome/Mus_musculus.GRCm39.106.gtf"
        " && rm genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
rule star_align:
    input:
        index = rules.star_index.output.dir,
        read1 = rules.fastp.output.read1,
        read2 = rules.fastp.output.read2
    output:
        bam = "starAlign/{sra}Aligned.sortedByCoord.out.bam",
        log = "starAlign/{sra}Log.final.out"
    log:
        out = "starAlign/{sra}_star.out",
        err = "starAlign/{sra}_star.err"
    params:
        prefix = "starAlign/{sra}"       
    threads:
        32
    shell:
        "STAR --runThreadN {threads} --genomeDir {input.index} --genomeLoad LoadAndKeep --readFilesIn {input.read1} {input.read2} --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outTmpDir /tmp/TMPDIR/{wildcards.sra} 1> {log.out} 2> {log.err} "

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
        "htseq-count {input.alignment} {input.gtf_1} -i gene_id -m union"

rule normalise_and_DE_analysis:
    input:
        rules.htseq_count.output.raw_count
    output:
        "AML_gene_lists.csv"
    script:
        "Analysis.Rmd"
        