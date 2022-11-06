#!/usr/bin/env bash

###################
# HELP is triggered with no arguments, or with -h or --help
if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ -z $1 ]
then
	echo "Argument 1: Input folder path containing basecalled FASTQ files"
	echo "Argument 2: Output folder, destination for clean FASTQ files"
	echo ""
	echo "Usage example:"
	echo "filter_quality_bbduk.sh 01_basecalled_demultiplexed/run001 02_bbduk/run001"
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

###################
# Clean reads using bbduk.sh
echo ""
echo ""
echo "Filtering reads for sequencing run "$INPUT" ..."
echo ""
cd $INPUT

for BC in barcode*.fastq.gz
do
	SAMPLE=${BC//.fastq.gz/}

	echo "-> Quality-filtering sample "$SAMPLE" ..."
	bbduk.sh -Xmx8g in=$BC out=$OUTPUT"/"$BC \
	qtrim=lr trimq=7 minlength=200 maq=10 qin=33 qout=33 ow \
	&>$OUTPUT"/"$SAMPLE".bbduk.log"
	echo "-> Quality-filtering complete for sample "$SAMPLE
	echo ""
	echo ""
done

echo "Quality filtering finished for sequencing run "$INPUT
echo ""
echo ""
echo ""
