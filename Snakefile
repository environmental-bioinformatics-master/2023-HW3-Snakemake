## DEFINE VARAIBLES 

SRR_FILES = glob_wildcards()
FASTQC = expand() # all of the FASTQC files you expect to be produced by your pipeline (forward and reverse)
TRIMMED = expand() # all of the trimmed files you expect to be produced by your pipeline (forward and reverse)
SALMON_INDEX = directory('<placeholder>') # name the directory where you expect your salmon index)
#SALMON_QUANT = expand() # all of the quantification files for each sample (ONE PER SRR, not one per forward and reverse
#SALMON_MERGE = '' # the path to the merged Salmon output


rule all:
    input: FASTQC, TRIMMED,


rule trim_galore_forward:
    input: ""
    output:
        fastqc = "",
        trimmed = ""
    conda:
        "envs/trim.yaml"
    shell:
        """
        mkdir -p outputs/trimqc_forward
        trim_galore -q 20 --phred33 --illumina --length 20 -stringency 3 --fastqc -o outputs/trimqc_forward {input}
        """
        
rule trim_galore_reverse:
    input: ""
    output:
        fastqc = "",
        trimmed = ""
    conda:
        "envs/trim.yaml"
    shell:
        """
        mkdir -p outputs/trimqc_reverse
        trim_galore -q 20 --phred33 --illumina --length 20 -stringency 3 --fastqc -o outputs/trimqc_reverse {input}
        """

#rule salmon_index:
#    input: "S_debilis_eye_assembly_clean.fasta"
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
