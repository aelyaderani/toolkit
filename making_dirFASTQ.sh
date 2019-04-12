#!/bin/bash
#SBATCH --job-name=CellRanger_Count
#SBATCH -N 1
#SBATCH -t 0-48:00
#SBATCH -n 1
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=aelyaderani@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/aelyaderani/Single_Cell_Analysis/Cell_Ranger/oeFiles/count_%j.out
#SBATCH -e /scratch/aelyaderani/Single_Cell_Analysis/Cell_Ranger/oeFiles/count_%j.err

declare projectN=X5SCR
declare flowcell=HG57WDSXX

#making cell ranger format directories (with correct layout)
mkdir /scratch/aelyaderani/temp_samples/$flowcell
mkdir /scratch/aelyaderani/temp_samples/Reports
mkdir /scratch/aelyaderani/temp_samples/Stats

#find a list of directories where fastq for the spicific flowcell are located.
find /liang/fastq_DoNotTouch/*$flowcell -name Sample*_$projectN_*_$flowcell_* > /scratch/aelyaderani/directorynamelist.txt
#take the text file, and repalce / with , to creat a csv
cat /scratch/aelyaderani/directorynamelist.txt | tr -s '/' ',' > /scratch/aelyaderani/directorynamelist.csv
#grabing only the 3rd column with the name needed for creating new directories
cut -d',' -f3 /scratch/aelyaderani/directorynamelist.csv > /scratch/aelyaderani/steponeDirectory.txt
cat /scratch/aelyaderani/steponeDirectory.txt | tr -s ‘_’ ',' > /scratch/aelyaderani/directorynamelist.csv
#cuting out the unnessory part of the name so cell ranger can recognize them.
cut -d',' -f2-10 /scratch/aelyaderani/directorynamelist.csv > /scratch/aelyaderani/steptwoDirectory.txt
#grabing LID list for the flowcell
cut -d',' -f9 /scratch/aelyaderani/directorynamelist.csv > /scratch/aelyaderani/LID.csv
cat /scratch/aelyaderani/steptwoDirectory.txt | tr -s ‘,’ ‘_’ > /scratch/aelyaderani/directorynamelist.csv

#removing duplicates from the csv list (LID and directory names)
sort -u /scratch/aelyaderani/directorynamelist.csv -o /scratch/aelyaderani/directorynamelist.csv
sort -u /scratch/aelyaderani/LID.csv -o /scratch/aelyaderani/LID.csv

#jump to flowcell directory we made earlier above and creat new directory from the directorynamelist.csv
cd /scratch/aelyaderani/temp_samples/$flowcell
mkdir $(</scratch/aelyaderani/directorynamelist.csv cut -d/ -f1)

#make variables with LID names
i=1
while IFS= read -r line; do
    declare primaryNode$((i++))="$line"
done </scratch/aelyaderani/LID.csv

count="$((i-1))"

#use the LID variables made to find and move the correct fastq files to the correct corresponding directory
for ((idx=1; idx<=count; idx++)); do
    temp=primaryNode$idx
    rsync -a --progress --stats --human-readable /liang/fastq_DoNotTouch/*$flowcell/*/*${!temp}*/*.fastq.gz  /scratch/aelyaderani/temp_samples/$flowcell/*${!temp}*
done
