SRA,FRR = glob_wildcards("rawReads/{sra}_{frr}.fastq.gz")

rule all:
    input:
        expand("rawQC/{sra}_{frr}_fastqc.{extension}", sra=SRA, frr=FRR,extension=["zip","html"]),
        expand("multiqc_report.html"),
        expand("trimmedreads{sra}_fastq.html", sra=SRA),
        
        
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
        rawqc="rawQC",
    output:
       "multiqc_report.html"
    shell:
        """
        multiqc rawQC
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
      
      
