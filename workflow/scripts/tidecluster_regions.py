#!/usr/bin/env python3
"""Read one TideCluster run (`TideCluster.py run_all` output directory,
tested with 1.21.3) into the tables the satellite cross-check uses, with
contig names translated to the pipeline's prepped-genome names.

Inputs, all `{dir}/{prefix}_*`:
  clustering.gff3        one feature per array, Name=TRC_n (required)
  tidehunter.gff3        every raw TideHunter array with its own monomer
  tidehunter_short.gff3  (consensus_length, copy_number); optional
  seqid_lengths.tsv      contig names/lengths of the assembly TideCluster ran
                         on; checked against the pipeline's fingerprint.tsv
                         (via name_map.tsv) so the coordinates are known to
                         refer to the same sequences (required)
  consensus_dimer_library.fasta   TRC consensi (required)
  tarean_report.tsv      per-TRC TAREAN monomer length, if a column for it
                         is recognisable (optional; logged)
  trc_superfamilies.csv  TRC -> superfamily (optional)
  rdna.tsv               TRCs TideCluster flags as rDNA (optional)

Writes:
  --regions  contig, start, end (0-based half-open), source (trc |
             tidehunter | tidehunter_short), id, monomer_len, copy_number
  --trc-info trc, n_arrays, array_bp, consensus_len, tarean_monomer_len,
             superfamily, rdna_flag
  --consensus  TRC consensi renamed <sample>:TRC_n
Parsing of the optional tables is deliberately tolerant (their layout
isn't documented): anything unrecognised is reported on stderr and left
as NA rather than guessed. Stdlib only."""

import argparse
import os
import re
import sys

from fasta_utils import iter_fasta, seq_id, write_fasta

TRC_RE = re.compile(r"\bTRC_\d+\b")
MISSING = {"", ".", "na", "nan", "none"}


def path_for(d, prefix, suffix, required=True):
    p = os.path.join(d, f"{prefix}_{suffix}")
    if not os.path.exists(p):
        if required:
            sys.exit(f"[tidecluster_regions] ERROR: missing TideCluster output {p}")
        print(f"[tidecluster_regions] note: optional {p} not found", file=sys.stderr)
        return None
    return p


def gff_attrs(field):
    out = {}
    for kv in field.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def read_gff(path):
    """Yield (seqid, start0, end, attrs) from a GFF3 (1-based inclusive)."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            yield f[0], int(f[3]) - 1, int(f[4]), gff_attrs(f[8])


def read_name_map(path):
    """TideCluster seqid (first token of the original header) -> prepped name."""
    m = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                m[seq_id(f[0])] = f[1]
    return m


def read_lengths(path, skip_hash=True):
    """name -> length from a 2-column TSV; headers/comment lines skipped."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if skip_hash and line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].strip().isdigit():
                out[f[0].strip()] = int(f[1])
    return out


def check_same_assembly(tc_lengths, name_map, fp_lengths):
    problems = []
    for tc_name, length in tc_lengths.items():
        ours = name_map.get(tc_name, tc_name)
        if ours not in fp_lengths:
            problems.append(f"{tc_name}: not in the pipeline genome")
        elif fp_lengths[ours] != length:
            problems.append(f"{tc_name}: length {length} vs pipeline {fp_lengths[ours]} ({ours})")
    mapped = {name_map.get(n, n) for n in tc_lengths}
    for ours in fp_lengths:
        if ours not in mapped:
            problems.append(f"{ours}: in the pipeline genome but not in TideCluster's seqid_lengths")
    return problems


def read_trc_table(path, label):
    """Rows of a small TSV/CSV that mention a TRC id: {TRC: [other fields]}.
    Returns ({}, None) if the file is absent."""
    rows = {}
    if not path:
        return rows
    with open(path) as fh:
        text = fh.read()
    sep = "," if path.endswith(".csv") else "\t"
    lines = [l for l in text.splitlines() if l.strip()]
    for line in lines:
        f = [x.strip() for x in line.split(sep)]
        ids = [x for x in f if TRC_RE.fullmatch(x)]
        if ids:
            rows[ids[0]] = [x for x in f if x != ids[0]]
    print(f"[tidecluster_regions] {label}: {len(lines)} lines, {len(rows)} with a TRC id"
          + (f"; first line: {lines[0]!r}" if lines else ""), file=sys.stderr)
    return rows


def read_tarean_monomers(path):
    """{TRC: monomer length} from the TAREAN report, if a column is recognisable."""
    if not path:
        return {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    print(f"[tidecluster_regions] tarean_report columns: {header}", file=sys.stderr)
    low = [h.lower() for h in header]
    i_len = next((i for i, h in enumerate(low) if "monomer" in h and ("len" in h or "size" in h)), None)
    if i_len is None:
        i_len = next((i for i, h in enumerate(low) if h in ("consensus_length", "monomer", "length")), None)
    if i_len is None:
        print("[tidecluster_regions] note: no monomer-length column recognised in the TAREAN "
              "report; tarean_monomer_len = NA", file=sys.stderr)
        return {}
    out = {}
    for r in rows:
        ids = [x for x in r if TRC_RE.search(x)]
        if ids and i_len < len(r):
            try:
                out[TRC_RE.search(ids[0]).group(0)] = int(float(r[i_len]))
            except ValueError:
                pass
    print(f"[tidecluster_regions] TAREAN monomer length from column '{header[i_len]}' "
          f"for {len(out)} TRCs", file=sys.stderr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--name-map", required=True)
    ap.add_argument("--fingerprint", required=True)
    ap.add_argument("--allow-seqid-mismatch", action="store_true")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--trc-info", required=True)
    ap.add_argument("--consensus", required=True)
    args = ap.parse_args()

    d, pre = args.dir, args.prefix
    name_map = read_name_map(args.name_map)
    tc_len = read_lengths(path_for(d, pre, "seqid_lengths.tsv"))
    fp_len = read_lengths(args.fingerprint)
    problems = check_same_assembly(tc_len, name_map, fp_len)
    if problems:
        msg = (f"TideCluster ran on a different assembly than the pipeline's {args.species} "
               f"({len(problems)} contig mismatches), e.g.:\n  " + "\n  ".join(problems[:10]))
        if not args.allow_seqid_mismatch:
            sys.exit("[tidecluster_regions] ERROR: " + msg +
                     "\nRerun TideCluster on the pipeline's prepped FASTA, or point the manifest at the "
                     "assembly TideCluster used. (satellite_crosscheck.allow_seqid_mismatch skips this "
                     "check; coordinates on mismatched contigs are then dropped.)")
        print("[tidecluster_regions] WARNING: " + msg, file=sys.stderr)
    print(f"[tidecluster_regions] {len(tc_len)} TideCluster contigs checked against the "
          f"pipeline genome: {len(problems)} mismatches", file=sys.stderr)

    def ours(seqid):
        """Pipeline contig name, or None if absent or a different length."""
        name = name_map.get(seqid, seqid)
        if name not in fp_len or tc_len.get(seqid, fp_len[name]) != fp_len[name]:
            return None
        return name

    regions = []
    trc_arrays, trc_bp = {}, {}
    dropped = 0
    for seqid, s, e, a in read_gff(path_for(d, pre, "clustering.gff3")):
        c = ours(seqid)
        trc = a.get("Name", "")
        if c is None or not TRC_RE.fullmatch(trc):
            dropped += 1
            continue
        regions.append((c, s, e, "trc", trc, "NA", "NA"))
        trc_arrays[trc] = trc_arrays.get(trc, 0) + 1
        trc_bp[trc] = trc_bp.get(trc, 0) + (e - s)
    for source, suffix in (("tidehunter", "tidehunter.gff3"), ("tidehunter_short", "tidehunter_short.gff3")):
        p = path_for(d, pre, suffix, required=False)
        if not p:
            continue
        for seqid, s, e, a in read_gff(p):
            c = ours(seqid)
            if c is None:
                dropped += 1
                continue
            regions.append((c, s, e, source, a.get("ID", "."),
                            a.get("consensus_length", "NA"), a.get("copy_number", "NA")))
    if dropped:
        print(f"[tidecluster_regions] dropped {dropped} features on unmatched contigs or without a TRC id",
              file=sys.stderr)

    with open(args.regions, "w") as out:
        out.write("contig\tstart\tend\tsource\tid\tmonomer_len\tcopy_number\n")
        for r in sorted(regions, key=lambda x: (x[0], x[1], x[2])):
            out.write("\t".join(map(str, r)) + "\n")

    cons_len = {}
    with open(args.consensus, "w") as out:
        for header, seq in iter_fasta(path_for(d, pre, "consensus_dimer_library.fasta")):
            m = TRC_RE.search(header)
            if not m:
                print(f"[tidecluster_regions] note: consensus without a TRC id skipped: {header[:60]}",
                      file=sys.stderr)
                continue
            trc = m.group(0)
            if trc in cons_len:
                continue
            cons_len[trc] = len(seq)
            write_fasta(out, f"{args.species}:{trc}", seq)

    monomer = read_tarean_monomers(path_for(d, pre, "tarean_report.tsv", required=False))
    superfam = read_trc_table(path_for(d, pre, "trc_superfamilies.csv", required=False), "trc_superfamilies")
    rdna = read_trc_table(path_for(d, pre, "rdna.tsv", required=False), "rdna")

    def trc_key(t):
        return int(t.split("_")[1])

    with open(args.trc_info, "w") as out:
        out.write("species\ttrc\tn_arrays\tarray_bp\tconsensus_len\ttarean_monomer_len\tsuperfamily\trdna_flag\n")
        for trc in sorted(set(trc_arrays) | set(cons_len), key=trc_key):
            sf = superfam.get(trc)
            out.write(f"{args.species}\t{trc}\t{trc_arrays.get(trc, 0)}\t{trc_bp.get(trc, 0)}\t"
                      f"{cons_len.get(trc, 'NA')}\t{monomer.get(trc, 'NA')}\t"
                      f"{sf[0] if sf else 'NA'}\t{'True' if trc in rdna else 'False'}\n")
    print(f"[tidecluster_regions] {len(trc_arrays)} TRCs, {sum(trc_arrays.values())} arrays, "
          f"{sum(trc_bp.values())} bp; {len(cons_len)} consensi", file=sys.stderr)


if __name__ == "__main__":
    main()
