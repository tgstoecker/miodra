configfile: "config.yaml"


rule all:
    input:
        "remapped_insertion_formatting.check"


# get coordinates -52/+51 of insertion and keep only insertions with at least two reads of support
# remove header line
# remove all lines that go into negative in second column
rule transform_insertions:
    input: 
        "/"
    output:
        "transformed_insertions.tsv"
    params:
        insertions=config["MuWU_insertion_file"],
    conda:
        "environment.yaml"
    shell:
        """
        awk -F "," '{{print $1,$2-51,$3+50}}' OFS="\\t" {params.insertions} | sort | uniq -d |
        tail -n +2 |
        awk 'BEGIN {{ OFS=FS="\\t" }} $2 !~ /^(-|0)/' > {output}
        """


# extract fasta sequences based on the prev. created file
rule bedtools_getfasta:
    input:
        "transformed_insertions.tsv"
    output:
        "insertions_up_downstream.tsv"
    params:
        genome=config["genome"],
    conda:
        "environment.yaml"
    shell:
        "bedtools getfasta -fi {params.genome} -bed {input} > {output}"


# align custom fasta reads to B73 reference
rule bowtie2_align:
    input:
        "insertions_up_downstream.tsv"
    output:
        "filtered.bam"
    params:
        index=config["bowtie2-index-basename"],
    threads:
        config["threads_bowtie_align"]
    conda:
        "environment.yaml"
    shell:
        """
        bowtie2 -p {threads} -f -x {params.index} -U {input} |\
        samtools view -b -F 4 > {output}
        """


# create from -bam file new file which contains:
# original coordinates in original mapping; chromosome of mapping, start+51; end+59 (so actual inserti$
# add header line
rule remapped_insertion_formatting:
    input:
        "filtered.bam"
    output:
        touch("remapped_insertion_formatting.check"),
    params:
        out_prefix=config["out_prefix"],
    conda:
        "environment.yaml"
    shell:
        """
        samtools view {input} |\
        awk -F "\\t" '{{print $1, $3, $4+51, $4+59}}' OFS="\\t" |\
        sed '1 i\\original_coords\\tchr\\tstart\\tend' > {params.out_prefix}.tsv
        """

