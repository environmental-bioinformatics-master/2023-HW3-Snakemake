# Snakemake Homework

In this homework, you will design and build a small snakemake workflow to parallel the homework on Transcriptomics. I have generated a smaller subset-ed version of the data used in Homework 3 in `/vortexfs1/omics/env-bio/collaboration/fastq_test` so that the data run more quickly-- but feel free to use the full data sets if you want. 

**The goal here is to:**
1. Trim the data using `trim_galore` 
2. Create an index of the reference transcriptome for `salmon`
2. Align the trimmed data against the reference transcriptome using `salmon quant`
3. Merge the counts produced in the 

First, you should see that you have a skeleton of a `Snakefile`. Take a look at what is already provided. At the end you should have a fully functional `Snakefile` that is able to recreate the first part of the analyses from the last homework with one command. 

This homework is designed to be done in order working command by command. Do not try to skip ahead! Work through one by one as each step builds on the previous. 

Take a look at the `Snakefile` that I provided for you.  As you can see I have sketched out the general structure and partially written two rules for you. You are going to be: 

1. Filling in the wildcards and expand statements to make `trim_galore` rule work. 
2. Writing a yaml file for salmon and fixing the shell command of the `salmon_index` rule. 
3. Writing (de novo) a rule for  `salmon_quant` and a rule for `salmon_merge`. 

## 0. Getting ready
First, open up tmux and spin up a srun session to work in (e.g. `srun -p scavenger --time=12:00:00 --ntasks-per-node=6 --mem=20gb --pty bash`)

Then, source activate your snakemake environment that you were using in class (likely `snakemake-class`). 

Finally, *link* (don't copy) all the data (`*fastq.gz` files) from the above directory into your local working environment in a folder called `raw_data`. So, your file structure should look like:

```bash
Homework6-Snakemake/
├── raw_data/
│   ├── SRR5004078_1.fastq.gz 
│   ├── SRR5004080_1.fastq.gz 
│   ├── SRR5004082_1.fastq.gz
│   ├── SRR5004083_1.fastq.gz 
│   ├── SRR5004085_1.fastq.gz 
│   ├── SRR5004090_1.fastq.gz 
│   ├── SRR5004091_1.fastq.gz 
│   ├── SRR5004093_1.fastq.gz 
│   ├── SRR5004094_1.fastq.gz 
│   ├── SRR5004095_1.fastq.gz 
│   ├── SRR5004096_1.fastq.gz 
│   └── SRR5004097_1.fastq.gz 
└── Snakefile
```

**Report the command you used to link the data:**
```

```

Now you are ready to start designing the pipeline. 

## 1. Trim Galore! 
Take a look at the `Snakefile` that was provided. You should see that you have one rule uncommented `rule trim_galore:`. To make this rule work you need to do a few things:
1. Define the `WILDCARDS` variable at the top to create a list of all the SRR files in the `raw_data/` directory.  
2. Use the `WILDCARDS`to design two `expand()` statements that create the output file lists for the fastqc output and trimmed fastq file from `trim_galore`. 
    **Note:** 
    - I want all files from this rule to be output into a directory called `outputs/trimqc/` so this should be in your expand statement.  
    - The tail of the fastqc output file looks like: `_trimmed_fastqc.zip`
    - The tail of the trimmed fastq file looks like: `_trimmed.fq.gz`
3. Fill in the appropriate `input` and `output` in the rule. 
4. Comment your rule please! 

It may be helpful to test out the wildcards and expand statements in a python interpreter like we did in class (with `from snakemake.io import *`). 

Once you have filled in the wildcards and expand statements you can start testing and refining this rule and the wildcards using the command: `snakemake --use-conda -p`. If it works, it should begin working through each of the files and creating the output. At the end of your run you should have a new set of files that look like this: 

```bash
outputs/
└── trimqc
    ├── SRR5004078_1.fastq.gz_trimming_report.txt
    ├── SRR5004078_1_trimmed_fastqc.html
    ├── SRR5004078_1_trimmed_fastqc.zip
    ├── SRR5004078_1_trimmed.fq.gz
    ├── SRR5004080_1.fastq.gz_trimming_report.txt
    ├── SRR5004080_1_trimmed_fastqc.html
    ├── SRR5004080_1_trimmed_fastqc.zip
    ├── SRR5004080_1_trimmed.fq.gz
    ├── SRR5004082_1.fastq.gz_trimming_report.txt
... etc.
```

Helpful hint: if you want to visualize file structure you can use the command `tree`. On the HPC it has to be loaded with `module load tree`. 

Copy and paste your final functioning `trim_galore` rule here:
```

```

> **QUESTION:** Here, many files are created by trim_galore (4 to be exact). What are the benefits of tracking all of them vs only some of them? Does it matter? 

> **ANSWER**: 


## 2. Salmon Index
 Once you have gotten the first rule working you are ready to move on to the `salmon` alignment and quantification!

The first step in `salmon` is creating an index that reads can be aligned again. This is done with the command:  
```bash
salmon index -t YOUR TRANSCRIPTOME -i SALMON_INDEX_DIRECTORY -k 25` 
```

Where `-t` indicates the path to your transcriptome fasta file and `-i` is the name of the index *directory* that will be created. Often, programs don't create a single file and it is easier to track a whole directory. For that you can use the call `directory()` command (see the value of `SALMON_INDEX`). 

For this rule, I have provided the appropriate input and output values. 

Here, you need to:
1. Uncomment the `salmon_index` rule. 
2. Fill in the shell command with appropriate wild cards. 
3. Create the yaml file (`envs/salmon.yaml`). Hint: take a look at the `fastqc.yaml` file to get an idea of how to structure it. 

Copy and paste your final functioning `salmon_index` rule here:
```

```

## 3. Salmon Quantification
Now, write a new rule called `salmon_quant` (there is a commented line in the current code that you can work off of). 

This rule should do the following: 
1. Run the command: 
```bash 
salmon quant -i SALMON_INDEX_DIRECTORY -l A -r SRR_TRIMMMED_FASTQ --validateMappings -o OUTPUT_DIRECTORY_SRR
```
2. Run the above command for *each* of the SRR samples. Each sample should be saved into a *directory* that is the SRR sample id. 
3.  It should depend on the creation of the salmon index file and the trimmed SRR data. 
4.  Each of the output directories should be located in `outputs/quant/`. 
5. Comment your rule! 

When it has worked, you should have a file structure that looks something like this:

```bash
├── quant
│   ├── SRR5004078
│   │   ├── aux_info
│   │   ├── cmd_info.json
│   │   ├── lib_format_counts.json
│   │   ├── libParams
│   │   ├── logs
│   │   └── quant.sf
│   ├── SRR5004080
│   ├── SRR5004082
│   ├── SRR5004083
│   ├── SRR5004085
│   ├── SRR5004090
│   ├── SRR5004091
│   ├── SRR5004093
│   ├── SRR5004094
│   ├── SRR5004095
│   ├── SRR5004096
│   └── SRR5004097
```
 Hints: It  might be helpful to look back at the `directory()` command above. Also note that the `directory()` command can be combined with  `expand()` command.  


Copy and paste your final functioning `salmon_quant` rule here:
```

```


## 4. Salmon Merge 
The command `quantmerge` in `salmon` can be used to combine your salmon quant tables from many different directories. Here, I want you to create a rule `salmon_merge` that creates a single file `outputs/quant/salmon_merged.counts`. 
1. The input should rely on *all* directories that were created in the Salmon Quant rule. 
2. The output should be one file: `outputs/quant/salmon_merged.counts`
3. It should run the command:
```
 salmon quantmerge --quants DIRECTORIES -o FINAL_QUANTFILE
```
One thing to consider here is that unlike previous rules (like `salmon_quant` and `trim_galore`) that wanted one sample at a time, this rule should take a list of all of the sample folders that you want to run this with. 

When it is done it should create a file that looks like this: 
```text
Name    SRR5004097      SRR5004094      SRR5004091      SRR5004083      SRR5004080      SRR5004085      SRR5004078      SRR5004090      SRR5004082      SRR5004095      SRR5004096      SRR5004093
contig#59353    0       0       0       0       0       0       0       0       0       0       0       0
contig#59351    0       0       0       0       0       0       0       0       0       0       0       43.9289
contig#59345    0       0       0       0       0       0       0       0       0       0       0       0
contig#59342    11.7669 0       0       0       0       0       0       0       0       0       0       0

```

Copy and paste your final functioning `salmon_merge` rule here:
```

```

##5. Snakemake-ing on `slurm`

Now, it is time to get `snakemake` to run `slurm` for you. To do this, I have provided a file called `cluster.yaml`. Open it up and take a look. This is the file that describes the resources required for running this file on slurm. The top `__default__: ` specifies default parameters. You don't need to change any of these except `account`. Here you should fill in your slurm username. But, you can see that some important info is provided (like run time in minutes, memory in GB). 

Below this you will see each of the rules you wrote specified. I have filled in one rule with requirements already for you (`trim_galore`). For the others, I leave it to you to fill in the requirements. Follow the formatting of the first rule, fill in the others. It might be helpful to look back at the suggested slurm requirements for the other rules in Homework 3.

Once you have saved the cluster file it is time to execute it on slurm. Generally, this is the command that you would use to do it (also saved in the file `submit_snakemake.sh`):

```bash
snakemake --jobs 20 --use-conda --cluster-config cluster.yaml --cluster "sbatch --parsable --partition={cluster.queue} --job-name=ENVBIO.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"
```

- `--jobs` specifies the max number of jobs to be running at once
- `--cluster-config` points to our cluster.yaml file and provides the parameters for the running of snakemake
- `--cluster` is a mad libs style fill in the blank that pulls information from the `cluster.yaml` file and generates the sbatch scripts

Run `rm -rf` on the `outputs/` folder to completely delete it. 

Now, make sure you are in `tmux` and are running the `snakemake` environment and copy and paste the above command into your terminal. This should start the processing of the snakemake jobs. 

Create a new pane with tmux (`ctrl+B "`) and type `squeue -u YOUR_USERNAME`. You should see many jobs running under your username. Save the output of squeue to a file called `slurmjobs.list`. 

Allow snakemake to run to fruition! When it has finished-- navigate into the hidden folder `.snakemake/log`. Here, you should see all the log files from everytime you ran snakemake. Copy the last snakemake log file to your main directory `Homework6-Snakemake` and call it `snakemake-final-run.log`. 

## You have finished! 
Congratulations! 

Please estimate how long this homework took you here: 

For the final submission, please add the following files to your git:
- `README.md`
- `Snakefile`
- `cluster.yaml`
- `snakemake-final-run.log`
- `slurmjobs.list`
- `output/quant/salmon_merged.counts`



