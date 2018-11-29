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
First thing is loading gatk and bwa that will sets the job options
Second, give the sample directory which is 'SAMPDIR' 
Then set up the variables:
    _R1 for read 1
    _R2 for read 2 
    _REF for reference
    _UBAM  for unaligned bam file
    _MET for unaligned.marked.txt
   _ABAM for aligned.bam file
After setting the variables and determing the path we go for the commandes 
This step has 3 commandes:
1. Command to Convert FASTQ reads to uBAM
2. command for MarkIlluminaAdapters:  Reads a SAM or BAM file and rewrites it with new adapter-trimming tags.
3. Command to Align with BWA


```bash
This is the command to launch the job:
sbatch -N 1 -p short -t 2880 --export SAMPDIR=P001/P001_N scripts/map_to_reference.sh : for normal tissu
sbatch -N 1 -p short -t 2880 --export SAMPDIR=P001/P001_T scripts/map_to_reference.sh : for tumor tissu

## Step 2: Mark duplicates

useful links: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146

https://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates

*the first thing to do to write the script is to find the requiered arguments by typing on the terminal: 
     _module load gatk/4.0.11.0
     _gatk MarkDuplicates --h
# used the requiered argument and then try to find the optional arguments needed by the help of the link above. 
This command was used to launch the markduplicates:
sbatch -N 1 -p short -t 2880 --export SAMPDIR=P001/P001_N scripts/mark_duplicates.sh : for normal tissu
sbatch -N 1 -p short -t 2880 --export SAMPDIR=P001/P001_T scripts/mark_duplicates.sh : for tumor tissu
```
## Step 3: Base recalibration

useful links: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
  https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
  https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
  
Tools involved: BaseRecalibrator, Apply Recalibration, AnalyzeCovariates (optional)

def: The base recalibration step is performed per-sample and consists of applying machine learning to detect and correct for patterns of systematic errors in the base quality scores, which are confidence scores emitted by the sequencer for each base. Base quality scores play an important role in weighing the evidence for or against possible variant alleles during the variant discovery process, so it's important to correct any systematic bias observed in the data.

*the first thing to do to write the script is to find the requiered arguments by typing on the terminal: 
     _module load gatk/4.0.11.0
     _gatk BaseRecalibrator --h
     
  *The requiered argument are:
 --input,-I:String             BAM/SAM/CRAM file containing reads  This argument must be specified at least once.
                              Required.

--known-sites:FeatureInput    One or more databases of known polymorphic sites used to exclude regions around known
                              polymorphisms from analysis.  This argument must be specified at least once. Required.

--output,-O:File              The output recalibration table file to create  Required.

--reference,-R:String         Reference sequence file  Required.


This command was used to launch the base recalibrator:
sbatch -N 1 -p short -t 2880 --export SAMPDIR=P001/P001_N scripts/base_recalibration.sh   : this is for normal tissu
sbatch -N 1 -p short -t 2880 --export SAMPDIR=P001/P001_T scripts/base_recalibration.sh   : this is for tumor tissu

