#!/usr/bin/env bash

###################
# HELP is triggered with no arguments, or with -h or --help
if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ -z $1 ]
then
	echo "Argument 1: Input folder path containing clean FASTQ files"
	echo "Argument 2: Output folder, destination for Canu de novo assemblies"
	echo "Argument 3: Number of CPU cores for minimap2"
	echo ""
	echo "Usage example:"
	echo "assemble_canu.sh 02_bbduk/run001 03_canu/run001 16"
	exit 0
fi

###################
# Exit if Input directory doesn't exist
if [ ! -d $1 ]
then
	echo "Input directory not found, verify the path "$1
	exit 0
fi

###################
# Create OUTPUT directory if it doesn't exist
mkdir -p $2

###################
# Set variables from command arguments
INPUT=(`realpath $1`)
OUTPUT=(`realpath $2`)
CORES=$3

###################
# Run Canu on basecalled and demultiplexed data, align reads to assembly after completion
echo ""
echo ""
echo "De novo assembly of samples in sequencing run "$INPUT" ..."
echo ""
cd $INPUT

# 	genomeSize=5k \ this was the original value

for BC in barcode*.fastq.gz
do
	SAMPLE=${BC//.fastq.gz/}

	# 1. Canu assembly
	echo "-> Assembling sample "$SAMPLE" ..."
	time canu \
	-d $OUTPUT"/"$SAMPLE"_canu" \
	-p $SAMPLE \
	-nanopore $BC \
	contigFilter="2 0 1.0 0.5 0" \
	corMaxEvidenceCoverageGlobal=20 \
	corMaxEvidenceCoverageLocal=20 \
	corMhapSensitivity=high \
	corMinCoverage=0 \
	corOutCoverage=600 \
	genomeSize=3k \
	maxInputCoverage=600 \
	minInputCoverage=4 \
	minOverlapLength=100 \
	minReadLength=200 \
	stopOnLowCoverage=4 \
	"batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" \
	&>$OUTPUT"/"$SAMPLE"_canu.log"
	echo "-> Assembly finished for sample "$SAMPLE

	# 2. Read mapping
	ASSEMBLY=$OUTPUT"/"$SAMPLE"_canu/"$SAMPLE".contigs.fasta"
	if [ -f $ASSEMBLY ]
	then
		echo "-> Mapping reads to assembly "$ASSEMBLY" ..."
		samtools faidx $ASSEMBLY
		minimap2 -t $CORES -x map-ont -a --MD --secondary=no \
		$ASSEMBLY $BC > $OUTPUT"/"$SAMPLE"_canu/"$SAMPLE".sam"
		samtools view -u $OUTPUT"/"$SAMPLE"_canu/"$SAMPLE".sam" | \
		samtools sort -@4 -o $OUTPUT"/"$SAMPLE"_canu/"$SAMPLE".bam"
		samtools index $OUTPUT"/"$SAMPLE"_canu/"$SAMPLE".bam"

		# 3. Assembly directory cleanup
		rm $OUTPUT"/"$SAMPLE"_canu/"$SAMPLE".sam"
		rm -rf $OUTPUT"/"$SAMPLE"_canu/"$SAMPLE".seqStore"
		rm -rf $OUTPUT"/"$SAMPLE"_canu/canu-scripts"
		rm -rf $OUTPUT"/"$SAMPLE"_canu/correction"
		rm -rf $OUTPUT"/"$SAMPLE"_canu/trimming"
		rm -rf $OUTPUT"/"$SAMPLE"_canu/unitigging"
		echo "-> Mapping completed for sample "$SAMPLE
	else
		rm -rf $OUTPUT"/"$SAMPLE"_canu"
		echo "-> Assembly FAILED, removing  directory: "$OUTPUT"/"$SAMPLE"_canu"
	fi
	echo ""
	echo ""
done

echo "De novo assembly completed for sequencing run "$INPUT
echo ""
echo ""
echo ""
