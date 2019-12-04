## DEFINE VARAIBLES 

WILDCARDS = glob_wildcards()
FASTQC = expand()
TRIMMED = expand()
SALMON_INDEX = directory('copepod_index')
#SALMON_QUANT = expand()
#SALMON_MERGE = ''


rule all:
    input: FASTQC, TRIMMED,


rule trim_galore:
    input: ""
    output:
        fastqc = "",
        trimmed = ""
    conda:
        "envs/trim.yaml"
    shell:
        """
        mkdir -p outputs/trimqc
        trim_galore -q 20 --phred33 --illumina --length 20 -stringency 3 --fastqc -o outputs/trimqc {input}
        """

#rule salmon_index:
#    input: "hw3_copepod_txm.fasta.gz"
#    output: SALMON_INDEX
#    conda:
#        "envs/salmon.yaml"
#    shell:
#        """
#        FILL IN 
#        """
#
#rule salmon_quant:
#
#
#
#rule salmon_merge:
#
#
#
#
