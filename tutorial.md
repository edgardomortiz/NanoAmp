# How to analyze Nanopore sequencing data with NanoAmp

<p align="right">
Version: v20250128<br>
Author: Gentaro Shigita
</p>

## 0. Preparation

Open the `Terminal` app and use the `mkdir` (= <u>m</u>a<u>k</u>e <u>dir</u>ectory) command to create a new directory (= folder) in which to run your analysis and store the results:

```sh
mkdir ~/Desktop/MyAnalysis
```

Use the `cd` (= <u>c</u>hange <u>d</u>irectory) command to move into the directory you just created:

```sh
cd ~/Desktop/MyAnalysis
```

Use the `pwd` (= <u>p</u>rint <u>w</u>orking <u>d</u>irectory) command to check that you are in the right place:

```sh
pwd
```

Use the `ln` (= <u>l</u>i<u>n</u>k) command to create a shortcut to the raw data in the current directory (raw data should be in the `minknow_data` directory on the `Desktop`):

```sh
ln -fns ../minknow_data/MySequencingRun 00_raw_data
```

Finally, activate the `nanoamp` environment so that you have access to all software needed for this pipeline.

```sh
conda activate nanoamp
```

> [!TIP]
> Did you perform live basecalling during the sequencing run? Then you can skip the next step.  
> However, you need to prepare your files in a specific directory structure (one FASTQ file per barcode) for the subsequent steps to work.
>
> To do so, run the following:
>
> ```sh
> mkdir 01_basecalled_demultiplexed
> for barcode in 00_raw_data/*/*/*/fastq_*/*; do
>   cat ${barcode}/*.fastq.gz >> 01_basecalled_demultiplexed/${barcode##*/}.fastq.gz
> done
> ```
>
> Proceed to `2. Quality filtering`

## 1. Basecalling & Demultiplexing

The first step is to convert the current signals acquired by the sequencer into DNA sequences (= reads) using `Guppy` (This process is called "basecalling"). At the same time, reads are sorted according to the barcode sequences attached during the library preparation (this process is called "demultiplexing"). These can be done with the `basecall_demultiplex_guppy.sh` script.

Check the usage:

```console
basecall_demultiplex_guppy.sh

Argument 1: Input folder path containing FAST5 files
Argument 2: Output folder, destination for FASTQ files
Argument 3: Nanopore barcoding kit used in library, (SQK-RPB004, SQK-PBK004, SQK-RBK110-96)

Usage example:
basecall_demultiplex_guppy.sh 00_raw_data/run001 01_basecalled_demultiplexed/run001 SQK-RBK110-96
```

Following the usage, your command should look like this:

```sh
basecall_demultiplex_guppy.sh 00_raw_data 01_basecalled_demultiplexed SQK-RBK110-96
```

This will take some time. You can check the progress in the `guppy_basecaller.log` file created in the output directory.

## 2. Quality filtering

As some of the reads are of poor quality and unreliable, the next step is to exclude them and retain decent-quality reads using [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/). This can be done with the `filter_quality_bbduk.sh` script.

Check the usage:

```console
filter_quality_bbduk.sh

Argument 1: Input folder path containing basecalled FASTQ files
Argument 2: Output folder, destination for clean FASTQ files

Usage example:
filter_quality_bbduk.sh 01_basecalled_demultiplexed/run001 02_bbduk/run001
```

Following the usage, your command should look like this:

```sh
filter_quality_bbduk.sh 01_basecalled_demultiplexed 02_bbduk
```

## 3. *De novo* assembly

The next step is to assemble the individual reads and reconstruct the amplicon(s) using [Canu](https://github.com/marbl/canu). This can be done with the `assemble_canu.sh` script.

Check the usage:

```console
assemble_canu.sh

Argument 1: Input folder path containing clean FASTQ files
Argument 2: Output folder, destination for Canu de novo assemblies
Argument 3: Expected total assembly length in kilobases (approximate sum of amplicon lengths)
Argument 4: Number of CPU cores for minimap2

Usage example:
assemble_canu.sh 02_bbduk/run001 03_canu/run001 3k 16
```

> [!IMPORTANT]
> `Canu` does not accept an expected assembly length (3rd argument) shorter than 1 kb. Put `1k` even if your amplicon is supposed to be shorter than 1 kb.

Following the usage, your command should look like this:

```sh
assemble_canu.sh 02_bbduk 03_canu 1k 16
```

## 4. Polish assembly

Since Nanopore sequencing is still less accurate than other sequencing technologies (e.g., Sanger, Illumina, and PacBio HiFi), particularly for homopolymers (e.g., AAAAA...), the next step is to infer and correct sequencing errors in the assembled contigs using [Medaka](https://github.com/nanoporetech/medaka) (This process is called "polishing"). This can be done with the `polish_raw_medaka.sh` script.

Check the usage:

```console
polish_raw_medaka.sh

Argument 1: Input folder path containing clean FASTQ files
Argument 2: Input folder path containing Canu de novo assemblies
Argument 3: Output folder, destination for polished assemblies

Usage example:
polish_raw_medaka.sh 02_bbduk/run001 03_canu/run001 04_medaka/run001
```

Following the usage, your command should look like this:

```sh
polish_raw_medaka.sh 02_bbduk 03_canu 04_medaka
```

## 5. Separate haplotypes (optional)

[Canu](https://github.com/marbl/canu) often merges two or more different but similar haplotypes into a single collapsed contig.
You can separate these collapsed haplotypes using [Longshot](https://github.com/pjedge/longshot). This may be of particular interest if you are working on hybrids and/or polyploids and can be done with the `haplotype_polished_longshot.sh` script.

Check the usage:

```console
haplotype_polished_longshot.sh

Argument 1: Input folder path containing clean FASTQ files
Argument 2: Input folder path containing polished assemblies
Argument 3: Output folder, destination for polished haploid assemblies
Argument 4: Number of CPU cores for minimap2

Usage example:
haplotype_polished_longshot.sh 02_bbduk/run001 04_medaka/run001 05_medaka_haploid 16
```

Following the usage, your command should look like this:

```sh
haplotype_polished_longshot.sh 02_bbduk 04_medaka 05_longshot 16
```

## 6. Identify samples (optional)

The `match_db.py` script can be used to assign a species identification to each haplotype. This script is also useful for renaming sequences according to the correspondence between barcode and sample name.

Check the usage:

```console
match_db.py

Argument 1: Input folder containing the polished haploid assemblies for the sequencing run
Argument 2: Output folder, destination for final annotated FASTA files and species identifications for the sequencing run
Argument 3: Path to barcodes file used in the sequencing run
Argument 4: Path to database file in FASTA format

Usage example:
match_db.py 05_medaka_haploid/run001 06_final_seqs_ids/run001 barcodes/run001.txt db/Cucurbitaceae.fasta
```

Following the usage, your command should look like this:

```sh
match_db.py 05_longshot 06_fastas barcodes.txt ../nanoamp_refs/Cucurbitaceae_5markers_v2.fasta
```
