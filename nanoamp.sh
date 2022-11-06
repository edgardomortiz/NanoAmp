RUN=$1
CORES=$2
SEQ_KIT="SQK-RBK110-96"
# SEQ_KIT="SQK-RBK004"
DB="references/Cucurbitaceae_5markers_v2.fasta"

basecall_demultiplex_guppy.sh $HOME/Desktop/minknow_data/$RUN 01_basecalled_data/$RUN $SEQ_KIT

filter_quality_bbduk.sh 01_basecalled_data/$RUN 02_bbduk/$RUN

assemble_canu.sh 02_bbduk/$RUN 03_canu/$RUN $CORES

polish_raw_medaka.sh 02_bbduk/$RUN 03_canu/$RUN 04_medaka/$RUN

haplotype_polished_longshot.sh 02_bbduk/$RUN 04_medaka/$RUN 05_medaka_haploid/$RUN $CORES

match_db.py 05_medaka_haploid/$RUN 06_final_seqs_ids_v2/$RUN barcodes/$RUN".txt" $DB
