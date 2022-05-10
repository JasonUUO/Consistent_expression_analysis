SRA,FRR = glob_wildcards("rawReads/{sra}_{frr}.fastq.gz")

rule all:
    input:
        expand("rawQC/{sra}_{frr}_fastqc.{extension}", sra=SRA, frr=FRR,extension=["zip","html"]),
        expand("multiqc_report.html"),
        expand("trimmedreads{sra}_fastq.html", sra=SRA),
        expand("genome/GRCh38.primary_assembly.genome.fa"),
        expand("genome/gencode.v29.annotation.gtf"),
        expand("starAlign/{sra}Aligned.sortedByCoord.out.bam", sra=SRA),
        
        
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
    "Download Genome sequence, primary assembly (GRCh38) "
    output:
        fa = "genome/GRCh38.primary_assembly.genome.fa"
    shell:
        "cd genome && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz && gunzip GRCh38.primary_assembly.genome.fa.gz"


rule get_genome_gtf:
    "Download gtf annotations corresponding to our genome"
    output:
        gtf = "genome/gencode.v29.annotation.gtf"
    shell:
        "cd genome && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz && gunzip gencode.v29.annotation.gtf.gz"


rule star_index:
    input:
        gtf = rules.get_genome_gtf.output.gtf,
        fa = rules.get_genome_fa.output.fa       
    output:
        dir = directory("starIndex"),
    threads:
        16
    shell:
        "mkdir {output.dir} && STAR  --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.dir} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} "

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
        16
    shell:
        "STAR --runThreadN {threads} --genomeDir {input.index} --genomeLoad LoadAndKeep --readFilesIn {input.read1} {input.read2} --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outTmpDir /tmp/TMPDIR/{wildcards.sra} 1> {log.out} 2> {log.err} "

