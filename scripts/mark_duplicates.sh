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

# Get the input aligned.bam
ABAM=${SAMPDIR}/aligned.bam
[[ ! -e "$ABAM" ]] && echo "$ABAM not found" && exit 1

# Set the name for the output files
MBAM=${SAMPDIR}/marked.bam
MET=${SAMPDIR}/marked.txt

echo "***INPUT***"
echo "Alignment:         ${ABAM}"
echo ""
echo "***OUTPUT***"
echo "Marked:            ${MBAM}"
echo "metrics:           ${MET}"
echo ""

##########################################################################################
# Step 1: Mark Duplicates
##########################################################################################
gatk MarkDuplicates \
    --INPUT ${ABAM} \
    --OUTPUT ${MBAM} \
    --METRICS_FILE ${MET} \
    --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900
