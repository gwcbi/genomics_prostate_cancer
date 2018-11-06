# Pilot Analysis

## Download References

The reference files we will use come from the [GATK Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).


```bash
mkdir -p refs
cd refs

# Download reference
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai

# Download BWT index
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa

# Download SNPs
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

cd ..
```

## Raw Data

Under the directory production a folder called P001 has been created refering to patient number 1
Under P001 2 folders have been created : P001_T for tumor and P001_N for normal 
Each folder (e.g P001_N) contains 2 files that end by R1.fastq.gz and R2.fastq.gz

## Step 1: Map to Reference

TODO: Describe what you are doing in this step

```bash
sbatch -N 1 -p short -t 2880 --export SAMPDIR=P001/P001_N scripts/map_to_reference.sh
```
