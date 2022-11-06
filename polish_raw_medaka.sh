#!/usr/bin/env bash

###################
# HELP is triggered with no arguments, or with -h or --help
if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ -z $1 ]
then
	echo "Argument 1: Input folder path containing clean FASTQ files"
	echo "Argument 2: Input folder path containing Canu de novo assemblies"
	echo "Argument 3: Output folder, destination for polished assemblies"
	echo ""
	echo "Usage example:"
	echo "polish_raw_medaka.sh 02_bbduk/run001 03_canu/run001 04_medaka/run001"
	exit 0
fi

###################
# Exit if Input directories don't exist
if [ ! -d $1 ]
then
	echo "Input directory not found, verify the path "$1
	exit 0
fi
if [ ! -d $2 ]
then
	echo "Input directory not found, verify the path "$2
	exit 0
fi

###################
# Create OUTPUT directory if it doesn't exist
mkdir -p $3

###################
# Set variables from command arguments
FASTQ_DIR=(`realpath $1`)
CANU_DIR=(`realpath $2`)
OUTPUT=(`realpath $3`)

###################
# Run Medaka on previously created BAM, or map reads first if BAM was not created
echo ""
echo ""
echo "Polishing de novo assemblies for sequencing run "$CANU_DIR" ..."
echo ""
cd $CANU_DIR

for ASM in *_canu
do
	SAMPLE=${ASM//_canu/}
	echo "-> Verifying sample "$SAMPLE" ..."
	ASSEMBLY=$CANU_DIR"/"$ASM"/"$SAMPLE".contigs.fasta"
	BAM=$CANU_DIR"/"$ASM"/"$SAMPLE".bam"
	if [ ! -f $BAM ]
	then
		FASTQ=$FASTQ_DIR"/"$SAMPLE".fastq.gz"
		echo "-> Mapping reads to assembly "$ASSEMBLY" ..."
		samtools faidx $ASSEMBLY
		minimap2 -x map-ont -a --MD --secondary=no \
		$ASSEMBLY $FASTQ > $CANU_DIR"/"$ASM"/"$SAMPLE".sam"
		samtools view -u $CANU_DIR"/"$ASM"/"$SAMPLE".sam" | \
		samtools sort -@4 -o $CANU_DIR"/"$ASM"/"$SAMPLE".bam"
		samtools index $CANU_DIR"/"$ASM"/"$SAMPLE".bam"
		rm $CANU_DIR"/"$ASM"/"$SAMPLE".sam"
		echo "-> Mapping completed for sample "$SAMPLE
	fi
	echo medaka consensus \
	$BAM $OUTPUT"/"$SAMPLE".hdf" \
	--model r941_min_sup_g507 \
	--threads 2 "&>"$OUTPUT"/"$SAMPLE"_medaka.log &&" \
	medaka stitch $OUTPUT"/"$SAMPLE".hdf" $ASSEMBLY $OUTPUT"/"$SAMPLE"_medaka.fasta" \
	"&>>"$OUTPUT"/"$SAMPLE"_medaka.log" \
	>> $OUTPUT"/polish_raw_medaka.params"
done

echo "-> Polishing assemblies in parallel ..."
parallel -j 12 -a $OUTPUT"/polish_raw_medaka.params"
rm $OUTPUT"/"*.hdf $OUTPUT"/"*.bed
echo ""
echo "Assembly polishing completed for sequencing run "$CANU_DIR" ..."
echo ""
echo ""
echo ""
