#! /bin/bash

umask 0002
module load gatk/4.0.11.0

export _JAVA_OPTIONS=-Xmx60g
[[ -e ./tmp ]] && export TMPDIR=./tmp

# Get SUBJECT from command line if provided
[[ -n $1 ]] && SUBJECT=$1

# Get the normal and tumor sample directories
NDIR=${SUBJECT}/${SUBJECT}_N
TDIR=${SUBJECT}/${SUBJECT}_T

# Get the normal and tumor BAM files
NBAM=${NDIR}/final.bam
[[ ! -e "$NBAM" ]] && echo "$NBAM not found" && exit 1

TBAM=${TDIR}/final.bam
[[ ! -e "$TBAM" ]] && echo "$TBAM not found" && exit 1

# Get the input Panel of Normals
PON="refs/1000g_pon.hg38.vcf.gz"
[[ ! -e "$PON" ]] && echo "$PON not found" && exit 1

# Get the input gnomad VCF (germline resource)
GVCF="refs/af-only-gnomad.hg38.vcf.gz"
[[ ! -e "$GVCF" ]] && echo "$GVCF not found" && exit 1

# Get the input reference FASTA
REF="refs/Homo_sapiens_assembly38.fasta"
[[ ! -e "$REF" ]] && echo "$REF not found" && exit 1

# Set the name for the output files
SVCF=${SUBJECT}/1_somatic_m2.vcf.gz
BOUT=${SUBJECT}/2_tumor_normal_m2.bam

echo "***INPUT***"
echo "Normal BAM:          ${NBAM}"
echo "Tumor BAM:           ${TBAM}"
echo "Panel of Normal:     ${PON}"
echo "Germline resource:   ${GVCF}"
echo "Reference FASTA:     ${REF}"
echo ""
echo "***OUTPUT***"
echo "Somatic VCF:         ${SVCF}"
echo "Output BAM:          ${BOUT}"
echo ""


##########################################################################################
# Step 1: Call somatic variants with MUTECT2
##########################################################################################
gatk Mutect2 \
    -R ${REF} \
    -I ${TBAM} \
    -tumor ${SUBJECT}_T \
    -I ${NBAM} \
    -normal ${SUBJECT}_N \
    -pon ${PON} \
    --germline-resource ${GVCF} \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    -O ${SVCF} \
    -bamout ${BOUT}
