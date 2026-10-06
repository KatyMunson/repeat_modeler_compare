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
  kite/monomer_size_top3_estimats.csv
                         KITE per-array founder period (tab-separated;
                         TRC_ID, seqid, start, end, array_length,
                         founder_period, ...). Unlike TideHunter's -P, which
                         limits which arrays are *found*, KITE re-measures each
                         found array's own period up to 10 kb (25 kb when
                         extended), so founders > -P are normal (optional)
  cmd_args.json          TideHunter's maximum repeat-unit length (-P, or
                         25000 with --long; TideCluster default 3000) ->
                         --params. Arrays of longer units are not in any TRC.

Writes:
  --regions  contig, start, end (0-based half-open), source (trc |
             tidehunter | tidehunter_short | kite), id, monomer_len, copy_number
             (kite rows: id = TRC, monomer_len = founder_period,
             copy_number = multiplicity)
  --trc-info trc, n_arrays, array_bp, consensus_len, tarean_monomer_len,
             kite_founder_median (array-length-weighted; the x-axis of
             TideCluster's cluster-overview plot), kite_founder_n,
             superfamily, rdna_flag
  --consensus  TRC consensi renamed <sample>:TRC_n
Parsing of the optional tables is deliberately tolerant (their layout
isn't documented): anything unrecognised is reported on stderr and left
as NA rather than guessed. Stdlib only."""

import argparse
import csv
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
    Parsed with csv: TideCluster's R tables quote fields that may hold
    newlines."""
    rows = {}
    if not path:
        return rows
    with open(path, newline="") as fh:
        records = [r for r in csv.reader(fh, delimiter="," if path.endswith(".csv") else "\t") if any(r)]
    for f in records:
        f = [x.strip() for x in f]
        ids = [x for x in f if TRC_RE.fullmatch(x)]
        if ids:
            rows[ids[0]] = [x for x in f if x != ids[0]]
    print(f"[tidecluster_regions] {label}: {len(records)} records, {len(rows)} with a TRC id"
          + (f"; first record: {records[0]!r}" if records else ""), file=sys.stderr)
    return rows


def read_tarean_monomers(path):
    """{TRC: monomer length} from the TAREAN report (tarean_report.R writes it
    with write.table(quote=TRUE); its Consensus column holds newlines inside
    the quotes, so it must be read with csv, not line by line)."""
    if not path:
        return {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = [h.strip() for h in next(reader, [])]
        rows = [r for r in reader if any(r)]
    print(f"[tidecluster_regions] tarean_report columns: {header}", file=sys.stderr)
    low = [h.lower() for h in header]
    i_len = next((i for i, h in enumerate(low) if "monomer" in h and ("len" in h or "size" in h)), None)
    if i_len is None:
        i_len = next((i for i, h in enumerate(low) if h in ("consensus_length", "monomer", "length")), None)
    i_trc = low.index("trc") if "trc" in low else None
    if i_len is None or i_trc is None:
        print("[tidecluster_regions] note: no TRC / monomer-length column recognised in the TAREAN "
              "report; tarean_monomer_len = NA", file=sys.stderr)
        return {}
    out = {}
    for r in rows:
        if max(i_len, i_trc) >= len(r) or not TRC_RE.fullmatch(r[i_trc].strip()):
            continue
        try:
            out[r[i_trc].strip()] = int(float(r[i_len]))
        except ValueError:
            pass
    print(f"[tidecluster_regions] TAREAN monomer length from column '{header[i_len]}' "
          f"for {len(out)} TRCs", file=sys.stderr)
    return out


def read_kite(path):
    """KITE per-array rows: [(seqid, start0, end, trc, founder, multiplicity, length)]."""
    rows = []
    if not path:
        return rows
    with open(path) as fh:
        header = [h.strip().strip('"') for h in fh.readline().rstrip("\n").split("\t")]
        need = ("TRC_ID", "seqid", "start", "end", "founder_period")
        if any(c not in header for c in need):
            print(f"[tidecluster_regions] note: KITE table lacks {[c for c in need if c not in header]}; "
                  f"columns are {header[:8]}...; not used", file=sys.stderr)
            return rows
        ix = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if len(f) < len(header):
                f += [""] * (len(header) - len(f))
            try:
                start, end = int(f[ix["start"]]), int(f[ix["end"]])
                founder = float(f[ix["founder_period"]])
            except ValueError:
                continue  # NA founder (below KITE's size threshold)
            mult = f[ix["multiplicity"]] if "multiplicity" in ix and f[ix["multiplicity"]] else "NA"
            length = end - start + 1
            rows.append((f[ix["seqid"]], start - 1, end, f[ix["TRC_ID"]], founder, mult, length))
    print(f"[tidecluster_regions] KITE: {len(rows)} arrays with a founder period", file=sys.stderr)
    return rows


def weighted_median(pairs):
    pairs = sorted(pairs)
    total = sum(w for _, w in pairs)
    acc = 0
    for v, w in pairs:
        acc += w
        if acc >= total / 2:
            return v
    return None


def tidehunter_max_period(path):
    """Longest repeat unit TideHunter looked for: TideCluster's -T/--tidehunter_arguments
    -P (default 3000), or 25000 with --long. Returns (period, how it was found)."""
    import json
    if not path:
        return 3000, "cmd_args.json missing; TideCluster default -P 3000 assumed"
    with open(path) as fh:
        try:
            args = json.load(fh)
        except ValueError:
            return 3000, "cmd_args.json unreadable; TideCluster default -P 3000 assumed"
    if not isinstance(args, dict):
        args = {}
    if args.get("long") is True:
        return 25000, "--long (three TideHunter rounds up to 25000)"
    th = args.get("tidehunter_arguments")
    m = re.search(r"-P\s+(\d+)", th) if isinstance(th, str) else None
    if m:
        return int(m.group(1)), f"tidehunter_arguments '{th}'"
    return 3000, "no -P in cmd_args.json; TideCluster default -P 3000 assumed"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--name-map", required=True)
    ap.add_argument("--fingerprint", required=True)
    ap.add_argument("--allow-seqid-mismatch", action="store_true")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--trc-info", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--params", required=True, help="sample, tidehunter_max_period, source")
    args = ap.parse_args()

    d, pre = args.dir, args.prefix
    name_map = read_name_map(args.name_map)
    tc_len = read_lengths(path_for(d, pre, "seqid_lengths.tsv"))
    fp_len = read_lengths(args.fingerprint)
    problems = check_same_assembly(tc_len, name_map, fp_len)
    if problems:
        msg = (f"TideCluster ran on a different assembly than the pipeline's {args.sample} "
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
    kite_by_trc = {}
    for seqid, s, e, trc, founder, mult, length in read_kite(
            path_for(d, pre, "kite/monomer_size_top3_estimats.csv", required=False)):
        c = ours(seqid)
        if c is None:
            dropped += 1
            continue
        regions.append((c, s, e, "kite", trc, f"{founder:g}", mult))
        kite_by_trc.setdefault(trc, []).append((founder, length))
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
            write_fasta(out, f"{args.sample}:{trc}", seq)

    monomer = read_tarean_monomers(path_for(d, pre, "tarean_report.tsv", required=False))
    superfam = read_trc_table(path_for(d, pre, "trc_superfamilies.csv", required=False), "trc_superfamilies")
    rdna = read_trc_table(path_for(d, pre, "rdna.tsv", required=False), "rdna")

    max_p, how = tidehunter_max_period(path_for(d, pre, "cmd_args.json", required=False))
    print(f"[tidecluster_regions] TideHunter max period {max_p} ({how})", file=sys.stderr)
    with open(args.params, "w") as out:
        out.write(f"sample\ttidehunter_max_period\tsource\n{args.sample}\t{max_p}\t{how}\n")

    for trc in list(monomer):
        cap = min(cons_len.get(trc, float("inf")), trc_bp.get(trc, float("inf")))
        if monomer[trc] > cap:
            print(f"[tidecluster_regions] WARNING: TAREAN monomer {monomer[trc]} for {trc} exceeds its "
                  f"consensus/array length ({cap:g}); set to NA", file=sys.stderr)
            del monomer[trc]

    def kmed(trc):
        m = weighted_median(kite_by_trc.get(trc, []))
        return "NA" if m is None else f"{m:g}"

    def trc_key(t):
        return int(t.split("_")[1])

    with open(args.trc_info, "w") as out:
        out.write("sample\ttrc\tn_arrays\tarray_bp\tconsensus_len\ttarean_monomer_len\t"
                  "kite_founder_median\tkite_founder_n\tsuperfamily\trdna_flag\n")
        for trc in sorted(set(trc_arrays) | set(cons_len), key=trc_key):
            sf = superfam.get(trc)
            out.write(f"{args.sample}\t{trc}\t{trc_arrays.get(trc, 0)}\t{trc_bp.get(trc, 0)}\t"
                      f"{cons_len.get(trc, 'NA')}\t{monomer.get(trc, 'NA')}\t"
                      f"{kmed(trc)}\t{len(kite_by_trc.get(trc, []))}\t"
                      f"{sf[0] if sf else 'NA'}\t{'True' if trc in rdna else 'False'}\n")
    print(f"[tidecluster_regions] {len(trc_arrays)} TRCs, {sum(trc_arrays.values())} arrays, "
          f"{sum(trc_bp.values())} bp; {len(cons_len)} consensi", file=sys.stderr)


if __name__ == "__main__":
    main()
