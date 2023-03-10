#!/usr/bin/env python3
base = "/mnt/disk1/PROJECTS/SURPRISE/eshi_i_myshi/READS/reads_180223/"
suffix_1 = "_1.fq.gz"
suffix_2 = "_2.fq.gz"

SAMPLES, = glob_wildcards(base +"{sample}"+suffix_1)

want_all = []
want_all.append(expand("results/kraken2/{sample}.report", sample=SAMPLES))
want_all.append(expand("results/hisat2/{sample}.sam", sample=SAMPLES))
want_all.append(expand("results/data/{sample}_R1_nohum.fastq", sample=SAMPLES))
want_all.append(expand("results/spades/{sample}.fasta", sample=SAMPLES))
want_all.append(expand("results/spades/filt_contigs/{sample}_l350.fasta", sample=SAMPLES))
want_all.append(expand("results/kraken2/contigs/{sample}.report", sample=SAMPLES))
want_all.append(expand("results/catbat/catbat_res_{sample}_c2c.txt", sample=SAMPLES))
want_all.append(expand("results/catbat/catbat_res_{sample}_true.names", sample=SAMPLES))
want_all.append(expand("results/kraken2/contigs/{sample}.contigs.taxo.txt", sample=SAMPLES))
want_all.append(expand("results/kraken2/contigs/{sample}.contigs_kraken2_anno.xlsx", sample=SAMPLES))
want_all.append(expand("results/catbat/catbat_res_{sample}.contigs.taxo.txt", sample=SAMPLES))
want_all.append(expand("results/catbat/catbat_res_{sample}.contigs_catbat_anno.xlsx", sample=SAMPLES))
want_all.append(expand("results/blast/contigs/{sample}_l350_blastres.out", sample=SAMPLES))
want_all.append(expand("results/blast/contigs/{sample}_l350_blastres_header_w_taxo.xlsx", sample=SAMPLES))
rule all:
    input: want_all

rule kraken2:
    input: read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2
    output: kraken2_report="results/kraken2/{sample}.report",
            kraken2_out="results/kraken2/{sample}.out"
    params: kraken2_db="/mnt/disk1/DATABASES/kraken2/pluspf"
    threads: 32
    shell: """
        kraken2 --threads {threads} --confidence 0.7 --db {params.kraken2_db} {input.read1} {input.read2} --use-names --report {output.kraken2_report} --output {output.kraken2_out}
    """
rule map_human:
    input: read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2
    output: "results/hisat2/{sample}.sam"
    params: ref="/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set"
    threads: 32
    shell: """
        hisat2 -x  {params.ref} -1 {input.read1} -2 {input.read2} -S {output} -p {threads}
        samtools flagstat {output} | awk -F "[(|%]" 'NR == 3 {{print $2}}' > {output}.stat
    """
rule filter_human:
    input: "results/hisat2/{sample}.sam"
    output: R1="results/data/{sample}_R1_nohum.fastq",
            R2="results/data/{sample}_R2_nohum.fastq"
    params: sam_unmapped="results/hisat2/{sample}.unmapped.sam",
    shell: """
        samtools view -b -f 12 -F 256 {input} > {params.sam_unmapped}
        bedtools bamtofastq -i {params.sam_unmapped} -fq {output.R1} -fq2 {output.R2}
       #rm {input}
       #rm {params.sam_unmapped}
    """
rule spades:
    input: sel_read1="results/data/{sample}_R1_nohum.fastq",
           sel_read2="results/data/{sample}_R2_nohum.fastq"
    output: res_file="results/spades/{sample}.fasta"
    params: outfolder="results/spades/{sample}"
    threads: 32
    shell: """
            spades.py --meta -t 16 -1 {input.sel_read1} -2 {input.sel_read2} -o {params.outfolder}
            cp {params.outfolder}/contigs.fasta {output.res_file}
            #mkdir -p all_spades
            #mv {params.outfolder} all_spades/
    """
rule contigs_filtering:
    input: contigs_raw="results/spades/{sample}.fasta"
    output: contigs_filt="results/spades/filt_contigs/{sample}_l350.fasta"
    shell: """
           seqkit seq -m 350 -g {input.contigs_raw} > {output.contigs_filt}
    """
rule contigs_kraken2:
    input: contigs_filt="results/spades/filt_contigs/{sample}_l350.fasta"
    output: kraken2_report="results/kraken2/contigs/{sample}.report",
            kraken2_out="results/kraken2/contigs/{sample}.out"
    params: kraken2_db="/mnt/disk1/DATABASES/kraken2/pluspf"
    threads: 32
    shell: """
        kraken2 --threads {threads} --confidence 0.7 --db {params.kraken2_db} {input.contigs_filt} --report {output.kraken2_report} --output {output.kraken2_out}
    """
rule catbat_contigs:
    input: contigs_filt="results/spades/filt_contigs/{sample}_l350.fasta"
    output: c2c= "results/catbat/catbat_res_{sample}_c2c.txt"
    params: ref_db="/mnt/disk1/DATABASES/CAT/CAT_prepare_20210107/2021-01-07_CAT_database",
            taxo_ref="/mnt/disk1/DATABASES/CAT/CAT_prepare_20210107/2021-01-07_taxonomy",
            prefix= "results/catbat/catbat_res_{sample}"
    threads: 32
    shell: """
          CAT contigs -c {input.contigs_filt} --index_chunks 1 --out_prefix {params.prefix} -n {threads} -d {params.ref_db} -t {params.taxo_ref} --force
          mv {params.prefix}.contig2classification.txt {output.c2c}
    """
rule catbat_human_names:
    input:c2c="results/catbat/catbat_res_{sample}_c2c.txt"
    output:true_names="results/catbat/catbat_res_{sample}_true.names"
    params: taxo_ref="/mnt/disk1/DATABASES/CAT/CAT_prepare_20210107/2021-01-07_taxonomy"
    shell: """
          CAT add_names -i {input.c2c} -o {output.true_names} -t {params.taxo_ref} --only_official
   """
rule kraken_contigs_2nd_taxo:
    input:  kraken2_contigs_out="results/kraken2/contigs/{sample}.out"
    output: kraken2_contigs_taxo="results/kraken2/contigs/{sample}.contigs.taxo.txt",
    shell: """
          grep 'C' {input.kraken2_contigs_out} | cut -f2,3 > {output.kraken2_contigs_taxo}
    """
rule kraken2_taxonomy_contigs:
    input: kraken2_spcs="results/kraken2/contigs/{sample}.contigs.taxo.txt"
    output: kraken2_contigs_taxo="results/kraken2/contigs/{sample}.contigs_kraken2_anno.xlsx"
    shell: """
          Rscript r_scripts/taxonomizr_kraken2.R {input.kraken2_spcs}
         # mkdir xlsx
         # mv {output.kraken2_contigs_taxo} xlsx/
    """
rule catbat_true_taxo:
   input:catbat_truenames="results/catbat/catbat_res_{sample}_true.names"
   output: catbat_taxo="results/catbat/catbat_res_{sample}.contigs.taxo.txt"
   shell: """
         cut -f1,4 {input.catbat_truenames} | awk '{{(n=split ($2,a,";"));print $1,"\t",a[n]}}' | sed 's/ \t /\t/g' | sed '1d' > {output.catbat_taxo}
   """
rule catbat_taxonomy_contigs:
    input: catbat_taxo="results/catbat/catbat_res_{sample}.contigs.taxo.txt"
    output: catbat_contigs_anno="results/catbat/catbat_res_{sample}.contigs_catbat_anno.xlsx"
    shell: """
          Rscript r_scripts/taxonomizr_catbat.R {input.catbat_taxo}
         # mkdir xlsx
         # mv {output.catbat_contigs_anno} xlsx/
    """
rule contigs_blast:
    input: contigs_filt="results/spades/filt_contigs/{sample}_l350.fasta"
    output: blast_out="results/blast/contigs/{sample}_l350_blastres.out",
            blast_out_wh="results/blast/contigs/{sample}_l350_blastres_header.out"
    params: ref="/mnt/disk1/blastdb/nt/nt"
    threads: 32
    shell: """
           blastn -query {input.contigs_filt} -db {params.ref} -max_target_seqs 1 -max_hsps 1 -num_threads {threads} -evalue 1e-5 -outfmt '6 qseqid sseqid length pident qcovs staxid stitle' -out {output.blast_out}
           cat blast_title.txt {output.blast_out} > {output.blast_out_wh}
    """
rule blast_make_pretty:
    input: blast_out_wh="results/blast/contigs/{sample}_l350_blastres_header.out"
    output: blast_pretty="results/blast/contigs/{sample}_l350_blastres_header_w_taxo.xlsx"
    shell: """
           Rscript r_scripts/blast_make_pretty.R {input.blast_out_wh}
          # mkdir xlsx
          # mv {output.blast_pretty} xlsx/
    """
