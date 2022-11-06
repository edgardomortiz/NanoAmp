#!/usr/bin/env bash

###################
# HELP is triggered with no arguments, or with -h or --help
if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ -z $1 ]
then
	echo "Argument 1: Input folder path containing clean FASTQ files"
	echo "Argument 2: Input folder path containing polished assemblies"
	echo "Argument 3: Output folder, destination for polished haploid assemblies"
	echo "Argument 4: Number of CPU cores for minimap2"
	echo ""
	echo "Usage example:"
	echo "haplotype_polished_longshot.sh 02_bbduk/run001 04_medaka/run001 05_medaka_haploid 16"
	exit 0
fi

###################
# Exit if Input directory or if reference doesn't exist
if [ ! -d $1 ]; then
	echo "Input directory not found, verify the path "$1
	exit 0
fi
if [ ! -d $2 ]; then
	echo "Input directory not found, verify the path "$2
	exit 0
fi

###################
# Create OUTPUT directory if it doesn't exist
mkdir -p $3

###################
# Set variables from command arguments
FASTQ_DIR=(`realpath $1`)
MEDAKA_DIR=(`realpath $2`)
OUTPUT=(`realpath $3`)
CORES=$4

###################
# Run longshot haplotyping on polished assembly
echo ""
echo ""
echo "Haplotyping polished de novo assemblies for sequencing run "$MEDAKA_DIR" ..."
echo ""
cd $MEDAKA_DIR
# 0. Remove previous results if they exist
rm *phased* *haploid* *longshot* *fai

for FASTA in *_medaka.fasta
do
	SAMPLE=${FASTA//_medaka.fasta/}
	samtools faidx $FASTA

	# 1. Read mapping
	echo "-> Mapping reads to "$FASTA" ..."
	minimap2 -t $CORES -x map-ont -a --MD --secondary=no \
	$FASTA $FASTQ_DIR"/"$SAMPLE".fastq.gz" \
	> $SAMPLE"_medaka.sam"
	samtools view -u $SAMPLE"_medaka.sam" | samtools sort -@4 -o $SAMPLE"_medaka.bam"
	samtools index $SAMPLE"_medaka.bam"
	rm $SAMPLE"_medaka.sam"
	echo "-> Mapping completed for sample "$SAMPLE

	# 2. Longshot haplotyping
	echo "-> Phasing sample "$SAMPLE" ..."
	longshot -b $SAMPLE"_medaka.bam" \
	-f $FASTA \
	-o $SAMPLE"_medaka_phased.vcf" \
	-O $SAMPLE"_medaka_phased.bam" \
	-c 20 \
	-C 50000 \
	-a 13 \
	-y 13 \
	-E 0.20 &>$SAMPLE"_medaka.longshot.log"
	if [ -f $SAMPLE"_medaka_phased.bam" ]
	then
		samtools index $SAMPLE"_medaka_phased.bam"
	fi
	rm $SAMPLE"_medaka.bam" $SAMPLE"_medaka.bam.bai"
	echo "-> Phasing completed for sample "$SAMPLE

	# 3. Extraction of phased contigs
	echo "-> Extracting phased contigs from sample "$SAMPLE" ..."
	bgzip $SAMPLE"_medaka_phased.vcf"
	tabix -p vcf $SAMPLE"_medaka_phased.vcf.gz"
	cat $FASTA | vcf-consensus $SAMPLE"_medaka_phased.vcf.gz" -s SAMPLE -H 1 > $SAMPLE"_medaka_phased.hap1.fasta"
	cat $FASTA | vcf-consensus $SAMPLE"_medaka_phased.vcf.gz" -s SAMPLE -H 2 > $SAMPLE"_medaka_phased.hap2.fasta"
	echo "-> Phased contigs extracted for sample "$SAMPLE

	# 4. Deduplicate homozygous contigs
	echo "-> Dereplicating phased assembly for sample "$SAMPLE" ..."
	ASSEMBLY=$SAMPLE"_medaka_haploid.fasta"
	sed 's/ .*//g' $SAMPLE"_medaka_phased.hap"?".fasta" | \
	dedupe.sh in=stdin.fasta uniquenames overwrite \
	out=$ASSEMBLY &>$SAMPLE"_medaka_phased.dedupe.log"
	rm $SAMPLE"_medaka_phased.hap1.fasta" $SAMPLE"_medaka_phased.hap2.fasta"
	echo "-> Phased assembly dereplicated "$ASSEMBLY

	# 5. Map reads to polished assembly
	echo "-> Mapping reads to haploid assembly "$ASSEMBLY" ..."
	minimap2 -t $CORES -x map-ont -a --MD --secondary=no \
	$ASSEMBLY $FASTQ_DIR"/"$SAMPLE".fastq.gz" \
	> $SAMPLE"_medaka_haploid.sam"
	BAM=$SAMPLE"_medaka_haploid.bam"
	samtools view -u $SAMPLE"_medaka_haploid.sam" | samtools sort -@4 -o $BAM
	samtools index $BAM
	rm $SAMPLE"_medaka_haploid.sam"
	echo "-> Mapping to haploid assembly completed for sample "$SAMPLE

	# 6. Run medaka on previous BAM
	echo medaka consensus \
	$BAM $SAMPLE".hdf" \
	--model r941_min_sup_g507 \
	--threads 2 "&>"$SAMPLE"_medaka_haploid.log &&" \
	medaka stitch $SAMPLE".hdf" $ASSEMBLY $OUTPUT"/"$SAMPLE".fasta" \
	"&>>"$SAMPLE"_medaka_haploid.log" \
	>> polish_haploid_medaka.params
	echo ""
	echo ""
done
echo "-> Polishing haploid assemblies in parallel ..."
parallel -j 12 -a polish_haploid_medaka.params
rm *.hdf $OUTPUT"/"*.bed
echo "-> Haploid assembly polishing completed for sequencing run "$MEDAKA_DIR" ..."

cd $OUTPUT
for FASTA in *.fasta
do
	SAMPLE=${FASTA//.fasta/}
	samtools faidx $FASTA
	# 7. Final mapping of reads to polished haploid assembly
	echo "-> Mapping reads to polished haploid assembly "$FASTA" ..."
	minimap2 -t $CORES -x map-ont -a --MD --secondary=no \
	$FASTA $FASTQ_DIR"/"$SAMPLE".fastq.gz" > $SAMPLE".sam"
	samtools view -u $SAMPLE".sam" | samtools sort -@4 -o $SAMPLE".bam"
	samtools index $SAMPLE".bam"
	echo "-> Mapping to haploid polished assembly completed for sample "$SAMPLE

	# 8. Caluclate coverage statistics for each contig
	echo "-> Computing coverage statistics for "$FASTA" ..."
	pileup.sh in=$SAMPLE".sam" out=$SAMPLE".covstats.txt" &>$SAMPLE".mapstats.txt"
	rm $SAMPLE".sam"
	echo "-> Coverage statistics computed for "$FASTA
	echo ""
	echo ""
done

echo "Haplotyping of polished assemblies completed for sequencing run "$MEDAKA_DIR" ..."
echo ""
echo ""
echo ""
