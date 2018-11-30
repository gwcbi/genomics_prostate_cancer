#! /bin/bash

umask 0002
module load gatk/4.0.11.0

export _JAVA_OPTIONS=-Xmx60g
[[ -e ./tmp ]] && export TMPDIR=./tmp

# Get SUBJECT from command line if provided
[[ -n $1 ]] && SUBJECT=$1

# Get the tumor sample directories	
TDIR=${SUBJECT}/${SUBJECT}_T

# Get the normal and tumor BAM files
TBAM=${TDIR}/final.bam
[[ ! -e "$TBAM" ]] && echo "$TBAM not found" && exit 1

# Get the input gnomad VCF (germline resource)
GVCF="refs/af-only-gnomad.hg38.vcf.gz"
[[ ! -e "$GVCF" ]] && echo "$GVCF not found" && exit 1

# Get the input somatic VCF
SVCF=${SUBJECT}/1_somatic_m2.vcf.gz
[[ ! -e "$SVCF" ]] && echo "$SVCF not found" && exit 1

# Set the name for the output files
PTABLE=${SUBJECT}/3_pileupsummaries.table
CTABLE=${SUBJECT}/4_calculatecontamination.table
FVCF=${SUBJECT}/5_filtered.vcf.gz

echo "***INPUT***"
echo "Tumor BAM:           ${TBAM}"
echo "gnomad VCF:          ${GVCF}"
echo "Somatic VCF:         ${SVCF}"
echo ""
echo "***OUTPUT***"
echo "Pileup table:        ${PTABLE}"
echo "Contamination table: ${CTABLE}"
echo "Filtered VCF:        ${FVCF}"
echo ""

##########################################################################################
# Step 1: Get pileup summaries
##########################################################################################
gatk GetPileupSummaries \
    -I ${TBAM} \
    -V ${GVCF} \
    -L ${GVCF} \
    -O ${PTABLE}

##########################################################################################
# Step 2: Calculate contamination
##########################################################################################
gatk CalculateContamination \
    -I ${PTABLE} \
    -O ${CTABLE}

##########################################################################################
# Step 3: Calculate contamination
##########################################################################################
gatk FilterMutectCalls \
    -L chr7:50626094-57291591  \
    -V ${SVCF} \
    --contamination-table ${CTABLE} \
    -O ${FVCF}
