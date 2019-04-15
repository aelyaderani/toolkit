#!/bin/bash
#SBATCH --job-name=CellRanger_Count
#SBATCH -N 1
#SBATCH -t 0-20:00
#SBATCH -n 1
#SBATCH --cpus-per-task 9
#SBATCH --mem=36G
#SBATCH --mail-user=aelyaderani@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/aelyaderani/temp_samples/oeFiles/count_%j.out
#SBATCH -e /scratch/aelyaderani/temp_samples/oeFiles/count_%j.err

module load cellranger/3.0.2
declare all_lanes=1
declare out=/scratch/aelyaderani/HG57WDSXX

declare cellranger_ref=/home/tgenref/homo_sapiens/grch38_hg38/tool_specific_resources/cellranger/refdata-cellranger-GRCh38-3.0.0

declare sample
declare cell_count=4000
declare chemistry=SC5P-PE
declare mem=34
declare cores=9

i=1
while IFS= read -r line; do
    declare primaryNode$((i++))="$line"
done </scratch/aelyaderani/directorynamelist.csv

count="$((i-1))"

for ((idx=1; idx<=count; idx++)); do

    temp=primaryNode$idx
    lib_id=${!temp}
    fastq_dir=/scratch/aelyaderani/temp_samples/HG57WDSXX/${!temp}

    cd $out
    echo -e "cellranger count --id=${lib_id}\n
                --transcriptome=${cellranger_ref}\n
                --fastqs=${fastq_dir}\n
                --sample=${sample}\n
                --expect-cells=${cell_count}\n
               	--chemistry=${chemistry}\n
                --localmem=${mem}\n
                --localcores=${cores}"

    cellranger count --id=${lib_id} \
                --transcriptome=${cellranger_ref} \
                --fastqs=${fastq_dir} \
                --sample=${sample} \
                --expect-cells=${cell_count} \
                --chemistry=${chemistry} \
                --localmem=${mem} \
                --localcores=${cores}

    if [ $? -eq 0 ]
    then
        echo -e "CellRanger Count: ${lib_id} PASSED"
    else
        echo -e "CellRanger Count: ${lib_id} FAILED"
    fi

done
