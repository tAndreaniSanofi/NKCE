#!/bin/sh
#SBATCH --job-name=BCR_Reconstruction      # job name
#SBATCH --array=1-8                        # Number of Job (according to the number of splitted files in the exported variables)
#SBATCH --nodes=2                          # nodes
#SBATCH -c 2                               # cores
#SBATCH --mem=6G                           # memory
#SBATCH --time=72:00:00                    # time
#SBATCH --error=velocyto.err                # error file name
#SBATCH --output=velocyto.out               # output file name
#SBATCH --mail-user=Tommaso.Andreani@sanofi.com  # email
#SBATCH --mail-type=ALL                          # type notification


### Bash script to be run on the workdir of magellan for an easy access of the data. Do not run the script on the DDS storage folders! it will take 10x more time to run each step.

#create the array to run in parallel the workflow
export fileName=`sed -n "$SLURM_ARRAY_TASK_ID"p sample_file_names.txt`

#activate the conda environment with velocyto
conda activate AIDA-ODS

#sort bam file by name
samtools sort -t CB -O BAM -@ 35 -o $fileName.cellsorted_possorted_genome_bam.bam $fileName.possorted_genome_bam.bam
samtools index  $fileName.cellsorted_possorted_genome_bam.bam;
velocyto run10x -m repeat_msk.gtf mypath_to_sorted_files/$fileName.cellsorted_possorted_genome_bam.bam /annotation/refdata-cellranger-mm10-1.2.0/genes/genes.gtf;

