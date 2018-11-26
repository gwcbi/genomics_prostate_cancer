#! /bin/bash

umask 0002
module load gatk/4.0.11.0

export _JAVA_OPTIONS=-Xmx60g
[[ -e ./tmp ]] && export TMPDIR=./tmp

# Get SAMPDIR from command line if provided
[[ -n $1 ]] && SAMPDIR=$1

# Extract SUBJECT and SAMPLE from SAMPDIR
SUBJECT=$(dirname $SAMPDIR)
SAMPLE=$(basename $SAMPDIR)

# Get the input marked.bam
MBAM=${SAMPDIR}/marked.bam
[[ ! -e "$MBAM" ]] && echo "$MBAM not found" && exit 1

# Get the input known VCF
KVCF="refs/Homo_sapiens_assembly38.dbsnp138.vcf"
[[ ! -e "$KVCF" ]] && echo "$KVCF not found" && exit 1

# Get the input reference FASTA
REF="refs/Homo_sapiens_assembly38.fasta"
[[ ! -e "$REF" ]] && echo "$REF not found" && exit 1

# Set the name for the output files
BQSROUT=${SAMPDIR}/recalibration_report.grp
ARBAM=${SAMPDIR}/final.bam

echo "***INPUT***"
echo "Alignment:           ${MBAM}"
echo "Known VCF:           ${KVCF}"
echo "Reference FASTA:     ${REF}"
echo ""
echo "***OUTPUT***"
echo "BQSR Table:          ${BQSROUT}"
echo "Analysis ready BAM:  ${ARBAM}"
echo ""


##########################################################################################
# Step 1: Base Quality Score Recalibration
##########################################################################################
gatk BaseRecalibrator \
    --reference ${REF} \
    --input ${MBAM} \
    --known-sites ${KVCF} \
    --output ${BQSROUT}

##########################################################################################
# Step 2: Apply BQSR
##########################################################################################
gatk ApplyBQSR \
    --reference ${REF} \
    --bqsr-recal-file ${BQSROUT} \
    --input ${MBAM} \
    --add-output-sam-program-record \
    --output ${ARBAM}
