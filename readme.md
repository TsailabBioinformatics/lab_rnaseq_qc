# Templates for RNA-Seq lab specific QC

## Instructions

### 1. Prepare genome reference

Use the `ref_prep.sh` to prepare index for genome reference. Make sure you change the genome's path. Then submit the job. It need a little bit more memory. 

### 2. Copy the fastq files

put all the fastq files in the `fastq` directory. Unzip the fastq files. Change the name into the format of `{sample_name}.R{1/2}.fastq`. Create a file list as shown in `file.list` containing `{sample_name}` only.

### 3. Submit the job

```bash
sbatch QC_only.sh
```

### 4. Get the report

Use get_trim_sum.py to aggregate the QC info.

```
python get_trim_sum.py ./clean/
```
