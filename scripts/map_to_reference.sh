#! /bin/bash

umask 0002
module load gatk/4.0.11.0
module load bwa/0.7.17

export _JAVA_OPTIONS=-Xmx60g
[[ -e /scratch ]] && TMPDIR=/scratch

# Get SAMPDIR from command line if provided
[[ -n $1 ]] && SAMPDIR=$1

# Extract SUBJECT and SAMPLE from SAMPDIR
SUBJECT=$(dirname $SAMPDIR)
SAMPLE=$(basename $SAMPDIR)

# Get the input R1 and R2
R1=$(ls ${SAMPDIR}/*R1.fastq.gz)
[[ ! -e "$R1" ]] && echo "$R1 not found" && exit 1
R2=$(ls ${SAMPDIR}/*R2.fastq.gz)
[[ ! -e "$R2" ]] && echo "$R2 not found" && exit 1

# Set the reference input name
REF=refs/Homo_sapiens_assembly38.fasta


# Set the name for output files
UBAM=${SAMPDIR}/unaligned.bam
MET=${SAMPDIR}/unaligned.marked.txt
ABAM=${SAMPDIR}/aligned.bam

echo "***INPUT***"
echo "Read 1:            ${R1}"
echo "Read 2:            ${R2}"
echo "Reference:         ${REF}"
echo ""
echo "***OUTPUT***"
echo "uBAM:              ${UBAM}"
echo "metrics:           ${MET}"
echo "aBAM:              ${ABAM}"
echo ""

##########################################################################################
# Step 1: Convert FASTQ reads to uBAM
##########################################################################################
gatk FastqToSam \
    -F1 ${R1} \
    -F2 ${R2} \
    -O ${SAMPDIR}/tmp.bam \
    -RG ${SAMPLE} \
    -SM ${SAMPLE} \
    -LB ${SUBJECT} \
    -PL illumina \
    --TMP_DIR $TMPDIR

##########################################################################################
# Step 2: MarkIlluminaAdapters 
##########################################################################################
gatk MarkIlluminaAdapters \
    -I ${SAMPDIR}/tmp.bam \
    -O ${UBAM} \
    -M ${MET}

##########################################################################################
# Step 3: Align with BWA
##########################################################################################
gatk SamToFastq \
    -I ${UBAM} \
    -F /dev/stdout \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \
    --INTERLEAVE true \
    --NON_PF true \
    --TMP_DIR $TMPDIR | \
bwa mem \
    -M \
    -t 7 \
    -p ${REF} \
    /dev/stdin | \
gatk MergeBamAlignment \
    --ALIGNED_BAM /dev/stdin \
    --UNMAPPED_BAM ${UBAM} \
    --OUTPUT ${ABAM} \
    -R ${REF} \
    --CREATE_INDEX true \
    --ADD_MATE_CIGAR true \
    --CLIP_ADAPTERS false \
    --CLIP_OVERLAPPING_READS true \
    --INCLUDE_SECONDARY_ALIGNMENTS true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ATTRIBUTES_TO_RETAIN XS \
    --TMP_DIR $TMPDIR


