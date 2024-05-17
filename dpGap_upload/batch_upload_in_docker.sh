#!/bin/bash

parent_directory="/home/data/real/dbGaP_Fetal/"

# Loop through each subdirectory (RNA,ATAC,CHIP,scRNA...)
for subfolder in "$parent_directory"/*; do
    if [ -d "$subfolder" ]; then
        # Loop through each .fastq.gz
        for file in "$subfolder"/*.fastq.gz; do
            if [ -e "$file" ]; then
                echo "$file"
		            /home/aspera/connect/bin/ascp -L /home/temp_log -i /home/aspera/connect/etc/aspera_tokenauth_id_rsa -Q -l 200m -k 1 $file asp-dbgap@gap-submit.ncbi.nlm.nih.gov:protected
            fi
        done
    fi
done
