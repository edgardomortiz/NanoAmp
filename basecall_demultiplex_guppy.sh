#!/usr/bin/env bash

###################
# HELP is triggered with no arguments, or with -h or --help
if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ -z $1 ]
then
	echo "Argument 1: Input folder path containing FAST5 files"
	echo "Argument 2: Output folder, destination for FASTQ files"
	echo "Argument 3: Nanopore barcoding kit used in library, (SQK-RPB004, SQK-PBK004, SQK-RBK110-96)"
	echo ""
	echo "Usage example:"
	echo "basecall_demultiplex_guppy.sh 00_raw_data/run001 01_basecalled_demultiplexed/run001 SQK-RBK110-96"
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
KIT=$3

###################
# Start basecalling
echo ""
echo ""
echo "Basecalling sequencing run "$INPUT" ..."
echo ""

# The next 2 parameters seem to be incompatible with --do_read_splitting and/or --overlap 100 (sup model)
# --detect_mid_strand_barcodes \
# --min_score_barcode_mid 60.0 \

# ############################
# # Run guppy on RTX3080
# echo "Basecalling sequencing run "$INPUT"..."
# time ~/software/ont-guppy_6.0.1_linux64/bin/guppy_basecaller \
# --device auto \
# --chunk_size 3000 \
# --chunks_per_runner 768 \
# --min_qscore 9 \
# --barcode_kits $KIT \
# --trim_barcodes \
# --trim_adapters \
# --trim_primers \
# --min_score_barcode_front 75.0 \
# --do_read_splitting \
# --calib_detect \
# --compress_fastq \
# --input_path $INPUT \
# --save_path $OUTPUT \
# --recursive \
# --config dna_r9.4.1_450bps_sup_plant.cfg &>$OUTPUT"/guppy_basecaller.log"

############################
# Run guppy on RTX2070 Super
# ~/software/ont-guppy_6.0.1_linux64/bin/guppy_basecaller \
# --trim_barcodes, barcode trimming is enabled by default in v6.3.8
# --min_qscore 9 \
# Latest guppy v6.3.8
time ~/software/ont-guppy/bin/guppy_basecaller \
--device auto \
--chunk_size 3000 \
--chunks_per_runner 72 \
--barcode_kits $KIT \
--trim_adapters \
--trim_primers \
--min_score_barcode_front 75.0 \
--do_read_splitting \
--calib_detect \
--compress_fastq \
--input_path $INPUT \
--save_path $OUTPUT \
--recursive \
--config dna_r9.4.1_450bps_sup_plant.cfg &>$OUTPUT"/guppy_basecaller.log"

###################
# Final file/folder cleanup
cd $OUTPUT"/pass"

for BC in *
do
	cd $BC
	cat *.fastq.gz > ../../$BC".fastq.gz"
	cd ..
done

cd ..
tar -czf fail.tar.gz fail
tar -czf guppy_basecaller_full_logs.tar.gz guppy_basecaller_log-*.log
rm -rf pass fail guppy_basecaller_log-*.log
echo "Basecalling finished for "$INPUT
echo ""
echo ""
echo ""
