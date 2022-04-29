SRA,FRR = glob_wildcards("rawReads/{sra}_{frr}.fastq.gz")

rule all:
    input:
        expand("rawQC/{sra}_{frr}_fastqc.{extension}", sra=SRA, frr=FRR, extension=["gz","html"]),
        expand("trimmedreads{sra}_fastq.html", sra=SRA),
        


rule rawFastqc:
    input:
        rawread="rawReads/{sra}_{frr}.fastq.gz",
    output:
        gz="rawQC/{sra}_{frr}_fastqc.gz",
        html="rawQC/{sra}_{frr}_fastqc.html",
    threads:
        1
    params:
        path="rawQC/",
    shell:
        """
        fastqc {input.rawread} --threads {threads} -o {params.path}
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
        6
    shell:
        """
         fastp --thread {threads} -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -h {output.report_html}
        """