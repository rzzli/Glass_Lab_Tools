# Steps to upload to dbGAP
### 1. build and run docker
#### - mount E:\fetal to /home/data
#### - use nobody as username
```bash
docker build -t dbgap:dev "E:\dd\dbgap\"
docker run -it --user nobody --mount type=bind,source="E:\fetal",target=/home/data dbgap:dev /bin/bash
```
### 2. download and install aspera, then configure token
```bash
#(in some dir)
wget https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_4.2.9.728_linux_x86_64.tar.gz
tar -zxvf ibm*
bash ibm*sh
export ASPERA_SCP_PASS=********************
```
### 3. run a single command to upload 
#### - make sure correct path for ascp
#### - make sure /home/temp_log is made for log
```bash
/home/aspera/connect/bin/ascp -L /home/temp_log -i /home/aspera/connect/etc/aspera_tokenauth_id_rsa -Q -l 200m -k 1 /home/data/temp/ATAC_F002_HMG086_Microglia.fastq.gz asp-dbgap@gap-submit.ncbi.nlm.nih.gov:test
```
