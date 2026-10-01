#!/usr/bin/env python3
"""Where does .out annotate bases that no .align alignment covers?
Diagnostic only, not run by the pipeline.
Usage: python3 workflow/scripts/diagnose_align_vs_out.py <x.fa.out> <x.fa.align>"""
import sys, collections

def merge(iv):
    iv.sort(); out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out

def uncovered(a, b):
    """bases of merged list a not in merged list b -> list of (s,e)"""
    res = []; j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] < cur: j += 1
        k = j
        while k < len(b) and b[k][0] <= e:
            if b[k][0] > cur: res.append((cur, b[k][0] - 1))
            cur = max(cur, b[k][1] + 1); k += 1
        if cur <= e: res.append((cur, e))
    return res

out_iv = collections.defaultdict(list); out_hits = collections.defaultdict(list)
for line in open(sys.argv[1]):
    f = line.split()
    if len(f) < 11 or not f[0].isdigit(): continue
    s, e = sorted((int(f[5]), int(f[6])))
    out_iv[f[4]].append((s, e)); out_hits[f[4]].append((s, e, f[10], f[9], line.rstrip()))

aln_iv = collections.defaultdict(list); n_hdr = 0; odd = collections.Counter(); odd_ex = []
for line in open(sys.argv[2]):
    f = line.split()
    if not f or not f[0].isdigit(): continue
    try:
        float(f[1]); s, e = sorted((int(f[5]), int(f[6])))
    except (ValueError, IndexError):
        odd["digit-led line, not a header"] += 1
        if len(odd_ex) < 3: odd_ex.append(line.rstrip())
        continue
    if not f[7].startswith("("):
        odd["header with field 8 not '(left)'"] += 1
        if len(odd_ex) < 3: odd_ex.append(line.rstrip())
    n_hdr += 1; aln_iv[f[4]].append((s, e))

print(f".out hits: {sum(map(len, out_iv.values()))}   .align headers: {n_hdr}")
for k, v in odd.items(): print(f"  {k}: {v}")
for x in odd_ex: print("   e.g.", x)

tot_out = tot_aln = tot_miss = 0
per_contig = []; miss_by_class = collections.Counter(); examples = []
for c, iv in out_iv.items():
    mo = merge(iv); ma = merge(aln_iv.get(c, []))
    bo = sum(e - s + 1 for s, e in mo); ba = sum(e - s + 1 for s, e in ma)
    miss = uncovered(mo, ma); bm = sum(e - s + 1 for s, e in miss)
    tot_out += bo; tot_aln += ba; tot_miss += bm
    per_contig.append((bm, c, bo, ba))
    if bm:
        for s, e, cls, rep, raw in out_hits[c]:
            for ms, me in miss:  # crude: attribute to hits overlapping missing segments
                if ms <= e and me >= s:
                    ov = min(e, me) - max(s, ms) + 1
                    miss_by_class[cls.split("/")[0]] += ov
                    if len(examples) < 8 and ov > 200: examples.append(raw)
print(f"\nunion bp  .out {tot_out:,}   .align {tot_aln:,}   .out-not-in-.align {tot_miss:,} ({100*tot_miss/tot_out:.1f}%)")
no_aln = [x for x in per_contig if x[3] == 0 and x[2] > 0]
print(f"contigs with .out hits but NO .align alignments: {len(no_aln)} carrying {sum(x[2] for x in no_aln):,} bp")
print("\ntop contigs by missing bp (missing, contig, out_bp, align_bp):")
for x in sorted(per_contig, reverse=True)[:10]: print("  ", x)
print("\nmissing bp by .out class (approx):")
for k, v in miss_by_class.most_common(): print(f"  {k:15s} {v:,}")
print("\nexample .out lines with no .align coverage:")
for x in examples: print("  ", x)
