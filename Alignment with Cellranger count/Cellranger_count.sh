#!/bin/bash
#SBATCH --job-name=cellranger          # Job name
#SBATCH --export=ALL                   # Exports all the results at once
#SBATCH --partition=medium             # Partition name
#SBATCH --time=3-00:00                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --signal=2                     # When a job is within sig_time seconds of its end time, send it the signal sig_num
#SBATCH --no-requeue                   # Specifies that the batch job should never be requeued under any circumstances
#SBATCH --cpus-per-task=16             # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=128G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=iclemente.1@alumni.unav.es # Direction email send

# Load CellRanger
module load CellRanger/6.1.1

# Set parameters
localcores=$SLURM_CPUS_PER_TASK
localmem=128


# Run CellRanger Count (Repeat this command per sample)
cellranger count \
    --id=run_count_p28_GFP \
    --sample=P28  \
    --fastqs=/home/iclemente.1/run_cellranger_count/p28_fastq  \
    --transcriptome=/home/iclemente.1/run_cellranger_count/Mus_musculus_genome_GFP_unmasked \
    --localcores="${localcores}" \
    --localmem="${localmem}" \
    --nosecondary                     # Will skip the automatic downstream clustering

# Run CellRanger Count (Repeat this command per sample)
cellranger count \
    --id=run_count_p31_GFP \
    --sample=P31  \
    --fastqs=/home/iclemente.1/run_cellranger_count/p31_fastq  \
    --transcriptome=/home/iclemente.1/run_cellranger_count/Mus_musculus_genome_GFP_unmasked \
    --localcores="${localcores}" \
    --localmem="${localmem}" \
    --nosecondary                       # Will skip the automatic downstream clustering