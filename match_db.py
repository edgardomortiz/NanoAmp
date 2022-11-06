#!/usr/bin/env python3
"""
Copyright 2022 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/NanoAmp

This file is part of NanoAmp. NanoAmp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. NanoAmp is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with NanoAmp. If
not, see <http://www.gnu.org/licenses/>.
"""

import gzip
import math
import multiprocessing
import subprocess
import sys
import time
import urllib
from pathlib import Path


# Reverse complement dictionary
REV_COMP_DICT = {
    "A": "T", "T": "A", "G": "C", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "r": "y", "y": "r", "s": "s", "w": "w",
    "K": "M", "M": "K", "B": "V", "V": "B", "k": "m", "m": "k", "b": "v", "v": "b",
    "D": "H", "H": "D", "N": "N", "d": "h", "h": "d", "n": "n",
    ".": ".", "-": "-", "?": "?"
}

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence
    """
    return "".join([REV_COMP_DICT[n] for n in seq[::-1]])


def fasta_to_dict(fasta_path):
    """
    Turns a FASTA file given with `fasta_path` into a dictionary. For example, for the sequence:
    ```text
    >k157_0 flag=1 multi=2.0000 len=274
    ATATTGATATTTCATAATAATAGTTTTTGAACTAAAAAGAAATTTTTCCTCCAATTATGTGGG
    ```
    Returns the dictionary:
    ```
    {'k157_0' : {
        'description': 'flag=1 multi=2.0000 len=274',
        'sequence': 'ATATTGATATTTCATAATAATAGTTTTTGAACTAAAAAGAAATTTTTCCTCCAAT...'
    }}
    ```
    """
    if f"{fasta_path}".endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    fasta_out = {}
    with opener(fasta_path, "rt") as fasta_in:
        seq  = ""
        name = ""
        desc = ""
        for line in fasta_in:
            line = line.strip("\n")
            if not line: continue
            if line.startswith(">"):
                if seq:
                    fasta_out[name] = {
                        "description": desc,
                        "sequence": seq,
                    }
                    seq = ""
                if len(line.split()) > 1:
                    name = line[1:].split()[0]
                    desc = " ".join(line.split()[1:])
                else:
                    name = line[1:].rstrip()
                    desc = ""
            else:
                seq += line
        if seq:
            fasta_out[name] = {
                "description": desc,
                "sequence": seq,
            }
    return fasta_out


def dict_to_fasta(
    in_fasta_dict, out_fasta_path, wrap=0, sort=False, append=False, write_if_empty=False
):
    """
    Saves a `in_fasta_dict` from function `fasta_to_dict()` as a FASTA file to `out_fasta_path`
    """
    if f"{out_fasta_path}".endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    if append is True:
        action = "at"
    else:
        action = "wt"
    if in_fasta_dict:
        if sort: in_fasta_dict = dict(sorted(in_fasta_dict.items(), key=lambda x: x[0]))
        with opener(out_fasta_path, action) as fasta_out:
            for name in in_fasta_dict:
                header = f'>{name} {in_fasta_dict[name]["description"]}'.strip()
                seq = in_fasta_dict[name]["sequence"]
                if wrap > 0:
                    seq_out = "\n".join([seq[i:i + wrap] for i in range(0, len(seq), wrap)])
                else:
                    seq_out = seq
                fasta_out.write(f'{header}\n{seq_out}\n')
    else:
        if write_if_empty:
            with opener(out_fasta_path, action) as fasta_out:
                fasta_out.write("")
    return out_fasta_path


def blat(blat_path, min_identity, min_coverage, size_tolerance, max_hits, target_path, query_path):
    """
    Extract matches of miscellaneous DNA sequences by comparing the assemblies to a set of
    references, these can be formatted as the proteins references '>sample-locus_name'
    """

    start = time.time()

    # Create output directory
    blat_dna_out_dir  = Path(target_path.parent)
    blat_dna_out_file = Path(blat_dna_out_dir, f"{target_path.stem}.blat.psl")
    blat_dna_log_file = Path(blat_dna_out_dir, f"{target_path.stem}.blat.log")
    id_report_tsv     = Path(blat_dna_out_dir, f"{target_path.stem}.sample_report.tsv")

    blat_cmd = [
        f"{blat_path}",
        "-t=dna",
        "-q=dna",
        "-noHead",
        f"-minIdentity={min_identity}",
        f"{target_path}",
        f"{query_path}",
        f"{blat_dna_out_file}"
    ]
    with open(blat_dna_log_file, "w") as blat_log:
        blat_log.write(f"NanoAmp' BLAT command:\n  {' '.join(blat_cmd)}\n\n")
    with open(blat_dna_log_file, "a") as blat_log:
        subprocess.run(blat_cmd, stdout=blat_log, stderr=blat_log)

    id_report = blat_misc_dna_psl_to_dict(blat_dna_out_file, min_coverage, size_tolerance,
                                          max_hits, target_path, query_path)
    with open(id_report_tsv, "wt") as tsv_out:
        for row in id_report:
            tsv_out.write("\t".join(row) + "\n")
    print(f'-> Identification report for {target_path.stem} saved to {id_report_tsv}')
    return


def blat_misc_dna_psl_to_dict(
        psl_path, min_coverage, size_tolerance, max_hits, target_path, query_path
    ):
    """
    Parse .psl from BLAT, assemble greedily the partial hits, and return the best set of hits if
    the reference contains more than a single sequence of the same type (analogous to the reference
    file Angiosperms353.FAA)
    """

    def calculate_psl_identity(
            q_end, q_start, t_end, t_start, q_inserts, matches, rep_matches, mismatches
    ):
        """
        Adapted for the case of DNA vs DNA search only from:
        https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl
        """
        millibad = 0
        q_ali_size = q_end - q_start
        t_ali_size = t_end - t_start
        ali_size = q_ali_size
        if t_ali_size < q_ali_size:
            ali_size = t_ali_size
        if ali_size <= 0:
            return millibad
        size_dif = abs(q_ali_size - t_ali_size)
        insert_factor = q_inserts
        total = matches + rep_matches + mismatches
        if total != 0:
            round_away_from_zero = 3 * math.log(1 + size_dif)
            if round_away_from_zero < 0:
                round_away_from_zero = int(round_away_from_zero - 0.5)
            else:
                round_away_from_zero = int(round_away_from_zero + 0.5)
            millibad = (1000 * (mismatches + insert_factor + round_away_from_zero)) / total

        return 100.0 - millibad * 0.1


    def determine_matching_region(
            q_size, q_start, q_end, t_size, t_start, t_end, t_strand, size_tolerance
    ):
        """
        Determine if a contig matches entirely the query, or if it is a partial hit. In cases of
        partial hits determine if it is a split hit (proximal, middle, or distal) or if the hit
        is partial and subsumed within a larger stretch of sequence unrelated to the query
        """
        region = ""
        if q_end - q_start >= q_size * (1 - size_tolerance):
            region = "full"
        elif (q_start <= q_size * size_tolerance
              and (q_size - q_end) <= q_size * size_tolerance):
            region = "full"
        elif (q_start >= q_size * size_tolerance
              and (q_size - q_end) >= q_size * size_tolerance):
            region = "middle"
        elif q_start <= q_size * size_tolerance:
            region = "proximal"
        elif (q_size - q_end) <= q_size * size_tolerance:
            region = "distal"
        else:  # hit is only partial and surrounded by a large proportion of unmatched sequence
            region = "wedged"
        # Now check if the flanks in the contig are too long to be part of a multi-part hit
        left_flank_too_long = t_start / (t_end - t_start) < size_tolerance
        right_flank_too_long = (t_size - t_end) / (t_end - t_start) < size_tolerance
        if region == "middle":
            if left_flank_too_long or right_flank_too_long:
                region = "wedged"
        elif ((region == "proximal" and t_strand == "+")
              or (region == "distal" and t_strand == "-")):
            if right_flank_too_long: region = "wedged"
        elif ((region == "proximal" and t_strand == "-")
              or (region == "distal" and t_strand == "+")):
            if left_flank_too_long: region = "wedged"

        return region


    def format_coords(starts, ends):
        """
        Given a list of starts [1,30,50] and a list of ends [12,38,100] produce a string for
        decorating the sequence description as '1-12,30-38,50-100', don't make them 1-based yet,
        that is done by the functions 'write_gff3' and 'write_fastas_and_report'
        """
        formatted_coords = []
        for s, e in zip(starts, ends):
            if s > e:
                formatted_coords.append(f'{e+1}-{s}')
            elif s == e:
                formatted_coords.append(f'{s}-{e}')
            else:
                formatted_coords.append(f'{s+1}-{e}')
        return ",".join(formatted_coords)


    dna_hits = {}
    hit = {}
    hit_num = 1
    with open(psl_path) as psl_in:
        for line in psl_in:
            cols = line.split()
            matches, mismatches, rep_matches = int(cols[0]), int(cols[1]), int(cols[2])
            q_inserts, t_strand, q_name, q_size = int(cols[4]), cols[8], cols[9], int(cols[10])
            t_base_inserts, t_name, t_size = int(cols[7]), cols[13], int(cols[14])
            q_start, q_end = [int(cols[11])], [int(cols[12])]
            t_start, t_end = [int(cols[15])], [int(cols[16])]

            coverage = (matches + rep_matches + mismatches) / q_size * 100
            identity = calculate_psl_identity(
                q_end[-1], q_start[0], t_end[-1], t_start[0],
                q_inserts, matches, rep_matches, mismatches
            )
            score = (matches + rep_matches - mismatches) / q_size
            wscore = score * ((matches + rep_matches + mismatches) / q_size)
            region = determine_matching_region(q_size, q_start[0], q_end[-1],
                                               t_size, t_start[0], t_end[-1],
                                               t_strand, size_tolerance)
            gapped = not bool(region == "full")

            # Ignore hits with not enough coverage of the reference, or partial hits immersed in
            # larger segments of unrelated sequence
            if not (coverage < min_coverage and region == "wedged"):

                # Compose hit record with columns from the .psl line
                # Many of these data won't make to final identification table, perhaps cleanup
                # unused fields later
                hit = {
                    "ref_name":     q_name,
                    "ref_taxon":    "-".join(q_name.split("-")[:-1]),
                    "ref_marker":   q_name.split("-")[-1],
                    "ref_size":     q_size,
                    "q_start":      q_start,                      # list of starts
                    "q_end":        q_end,                          # list of ends
                    "match_len":    matches + rep_matches + mismatches,
                    "hit_id":       f"DNA{hit_num}",
                    "hit_contig":   t_name,
                    "t_start":      t_start,                      # list of starts
                    "t_end":        t_end,                          # list of ends
                    "match_coords": format_coords(t_start, t_end),
                    "strand":       t_strand,
                    "matches":      matches + rep_matches,
                    "mismatches":   mismatches,
                    "coverage":     coverage,
                    "identity":     identity,
                    "score":        score,
                    "wscore":       wscore,
                    "region":       region,
                    "gapped":       gapped,
                }
                hit_num += 1
                if t_name not in dna_hits:
                    dna_hits[t_name] = [dict(hit)]
                else:
                    dna_hits[t_name].append(dict(hit))

    # Sort matches to each contig from highest to lowest 'score'
    for contig in dna_hits:
        dna_hits[contig] = sorted(dna_hits[contig], key=lambda i: i["score"], reverse=True)

    # Load target (assembly) and extract depth and voucher data
    target_fasta = fasta_to_dict(target_path)
    target_data = {}
    for contig in target_fasta:
        contig_data = {}
        for info in target_fasta[contig]["description"].replace("[", "").replace("]", "").split():
            contig_data[info.split("=")[0]] = info.split("=")[1]
        target_data[contig] = contig_data

    # Load query (reference), assuming every sequence includes marker name after last "-"
    query_fasta = fasta_to_dict(query_path)
    marker_names = sorted(list(set([seq_name.split("-")[-1] for seq_name in query_fasta])))

    # Collate final identification table
    id_report = [[
        "voucher",
        "run",
        "barcode",
        "marker",
        "contig_name",
        "contig_length",
        "depth",
        "depth_sd",
        "ranking",
        "ref_taxon",
        "score",
        "identity",
        "coverage",
        "match_coords",
        "strand"
    ]]
    for marker in marker_names:
        found=False
        for contig in dna_hits:
            for i in range(len(dna_hits[contig])):
                if marker == dna_hits[contig][i]["ref_marker"] and i < max_hits:
                    found=True
                    id_report.append([
                        f'{target_data[contig]["voucher"]}',
                        f'{target_data[contig]["run"]}',
                        f'barcode{target_data[contig]["barcode"]}',
                        f'{marker}',
                        f'{contig}',
                        f'{target_data[contig]["length"]}',
                        f'{target_data[contig]["depth"]}',
                        f'{target_data[contig]["depth_sd"]}',
                        f'{i+1}',
                        f'{dna_hits[contig][i]["ref_taxon"]}',
                        f'{dna_hits[contig][i]["score"]:.4f}',
                        f'{dna_hits[contig][i]["identity"]:.2f}',
                        f'{dna_hits[contig][i]["coverage"]:.2f}',
                        f'{dna_hits[contig][i]["match_coords"]}',
                        f'{dna_hits[contig][i]["strand"]}',
                    ])
        if not found:
            id_report.append([
                f'{target_data[contig]["voucher"]}',
                f'{target_data[contig]["run"]}',
                f'barcode{target_data[contig]["barcode"]}',
                f'{marker}',
                f'NA',
                f'NA',
                f'NA',
                f'NA',
                f'NA',
                f'NA',
                f'NA',
                f'NA',
                f'NA',
                f'NA',
                f'NA',
            ])

    return id_report


if __name__ == '__main__':

    # Quick program help
    if len(sys.argv) < 5:
        print("Argument 1: Input folder containing the polished haploid assemblies for the sequencing run")
        print("Argument 2: Output folder, destination for final annotated FASTA files and species identifications for the sequencing run")
        print("Argument 3: Path to barcodes file used in the sequencing run")
        print("Argument 4: Path to database file in FASTA format")
        print()
        print("Usage example:")
        print("match_db.py 05_medaka_haploid/run001 06_final_seqs_ids/run001 barcodes/run001.txt db/Cucurbitaceae.fasta")
        sys.exit(0)

    # Parse input arguments
    input_dir     = Path(sys.argv[1])
    output_dir    = Path(sys.argv[2])
    barcodes_file = Path(sys.argv[3])
    fasta_db_file = Path(sys.argv[4])

    # Define constants, hard-coded arguments, changed inside script
    blat_path      = "blat" # change to specific BLAT to path if needed
    min_identity   = 75     # minimum identity to report matches
    min_coverage   = 20     # minimum recovered percentage to report matches
    size_tolerance = 0.05   # proportion of match size tolerated to still consider full match
    max_hits       = 5      # number of top hits to report in full individual sample reports
    cpus           = 'all'  # change to 'one' for debugging

    # Verify input files and directories
    if not input_dir.exists():
        print(f'Input directory not found, please verify the given path {input_dir}')
        sys.exit(0)
    else:
        input_dir = input_dir.resolve()

    if not barcodes_file.exists():
        print(f'Barcodes file not found, please verify the given path {barcodes_file}')
        sys.exit(0)
    else:
        barcodes_file = barcodes_file.resolve()

    if not fasta_db_file.exists():
        print(f'FASTA Database file not found, please verify the given path {fasta_db_file}')
        sys.exit(0)
    else:
        fasta_db_file = fasta_db_file.resolve()

    if output_dir.exists():
        print(f'Output directory already exists, delete or rename {output_dir} to avoid overwriting previous results')
        sys.exit(0)
    else:
        output_dir.mkdir(parents=True)
        output_dir = output_dir.resolve()

    # List FASTA assemblies in input directory
    fastas = list(input_dir.glob("*.fasta"))

    # Create barcodes dictionary
    barcodes = {}
    with open(barcodes_file, "rt") as bc:
        for line in bc:
            line = line.strip("\n").split()
            if not line:
                continue
            else:
                barcodes[line[0]] = "_".join(line[1:])

    # Print arguments
    print()
    print(f'Input directory:     {input_dir}')
    print(f'Output directory:    {output_dir}')
    print(f'Barcodes file:       {barcodes_file}')
    print(f'FASTA database file: {fasta_db_file}')
    print()
    print(f'Found {len(fastas)} FASTA assemblies and {len(barcodes)} barcodes')
    print()
    print()

    # Rename FASTA assemblies, rename headers, include useful data in descriptions
    print(f'STEP 1: Renaming and adding sequencing depth info to FASTA assemblies')
    print()
    for fasta in fastas:
        barcode = fasta.stem
        covfile = Path(fasta.parent, f'{barcode}.covstats.txt')
        covstats = {}
        with open(covfile, "rt") as stats:
            for line in stats:
                if not line.startswith("#"):
                    s = line.strip("/n").split()
                    covstats[s[0]] = f'[depth={float(s[1]):.1f}] [depth_sd={float(s[10]):.1f}] [length={s[2]}]'
        seq_run = fasta.parts[-2]
        bc = barcode.replace("barcode", "BC")
        fasta_in = fasta_to_dict(fasta)
        fasta_out = {}
        for seq_name in sorted(fasta_in):
            if seq_name.endswith("_dd0"):
                hap = "B"
            elif seq_name.endswith("_dd1"):
                print(seq_name)
            else:
                hap = "A"
            new_seq_name = f'{seq_name.replace("000000", "").replace("_dd0", "")}{hap}'
            description = f'{covstats[seq_name]} [run={seq_run}] [barcode={bc.replace("BC", "")}] [voucher={barcodes[barcode]}]'
            fasta_out[new_seq_name] = {
                "sequence": fasta_in[seq_name]["sequence"],
                "description": description
            }
        fastaout_path = Path(output_dir, f'{barcodes[barcode]}__{bc}__{seq_run}.fasta')
        dict_to_fasta(fasta_out, fastaout_path, wrap=100)
        print(f'-> Final FASTA output for {fasta.stem} saved to {fastaout_path}')


    # Prepare list of parameters for each assembly for running BLAT
    fastas = list(output_dir.glob("*.fasta"))
    blat_params = []
    for fasta in fastas:
        blat_params.append((
            blat_path,
            min_identity,
            min_coverage,
            size_tolerance,
            max_hits,
            fasta,         # target: the haploid assembly
            fasta_db_file, # query : the FASTA database to search
        ))

    # Run in BLAT in parallel if using multiple cores, or serial if using just 1
    print()
    print()
    print(f'STEP 2: Comparing assembled sequences to refence database')
    print()
    if cpus == "all":
        process = multiprocessing.Pool(multiprocessing.cpu_count())
        for i in range(len(fastas)):
            process.apply_async(blat, blat_params[i])
        process.close()
        process.join()
    else:
        for i in range(len(fastas)):
            blat(*blat_params[i])

    # Collect all results for sequencing run in a single table, only keep top hit for each contig
    print()
    print()
    print(f'STEP 3: Sequencing run identification report')
    print()
    sample_reports = sorted(list(output_dir.glob("*.sample_report.tsv")))
    run_report = [[
        "voucher",
        "run",
        "barcode",
        "marker",
        "contig_name",
        "contig_length",
        "depth",
        "depth_sd",
        "ranking",
        "ref_taxon",
        "score",
        "identity",
        "coverage",
        "match_coords",
        "strand"
    ]]
    for sample_report in sample_reports:
        with open(sample_report, "rt") as rin:
            fasta_in = fasta_to_dict(Path(f'{sample_report}'.replace(".sample_report.tsv", ".fasta")))
            fasta_out = {}
            for line in rin:
                if not line.startswith("voucher"):
                    record = line.strip().split()
                    if record[8]:# in ["NA", "1"]:
                        run_report.append(record)
                    if record[8] != "NA":
                        new_seq_name = f'{record[0]}_{record[2].replace("barcode", "bc")}{record[4]}-{record[3]}'
                        if record[14] == "-":
                            seq = reverse_complement(fasta_in[record[4]]["sequence"])
                        else:
                            seq = fasta_in[record[4]]["sequence"]
                        fasta_out[new_seq_name] = {
                            "sequence": seq,
                            "description": fasta_in[record[4]]["description"],
                        }
            dict_to_fasta(fasta_out, Path(f'{sample_report}'.replace(".sample_report.tsv", ".fasta")))
    run_report_tsv = Path(output_dir, f'{output_dir.parts[-1]}.run_report.tsv')
    with open(run_report_tsv, "wt") as rout:
        for row in run_report:
            rout.write("\t".join(row) + "\n")
    print(f"Sequencing run identification report saved to {run_report_tsv}")
    print()
    print("Done!\n")