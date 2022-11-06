# NanoAmp
Species identification via Nanopore sequencing of barcoded amplicons


# Installation:
NanoAmp was tested on Linux, but should work on Mac as well. The dependencies can be installed in a separate conda environment with the following command:
```bash
conda create -n nanoamp -c bioconda -c conda-forge canu=2.2 parallel medaka=1.6.0 longshot tabix perl-vcftools-vcf bbmap ucsc-blat
conda activate nanoamp
```

Afterwards, this repository must be cloned:
```bash
git clone https://github.com/edgardomortiz/NanoAmp
```

Finally, add the directory where you cloned NanoAmp to your system $PATH


# Usage:
All steps can be run with a single command:
```bash
nanoamp.sh RUN_NAME CORES
```
Where `RUN_NAME` is the name of the sequencing run, assigned in MinKNOW. `CORES` indicates the number of processors to use.  
  
However, if you wan to run steps manually, the order of execution is:

1. Basecalling: `basecall_demultiplex_guppy.sh`
2. Quality filtering: `filter_quality_bbduk.sh`
3. Assembly: `assemble_canu.sh`
4. Polish assembly: `polish_raw_medaka.sh`
5. Separate haplotypes: `haplotype_polished_longshot.sh`
6. Identify samples: `match_db.py`

For help on each script just execute the script without extra arguments.
