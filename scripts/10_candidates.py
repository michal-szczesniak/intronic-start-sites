#!/usr/bin/env python3
"""
iSTART, Stage 0, step 1: candidate sense iRNAs + control sets.

Input : GENCODE GTF (e.g. gencode.v47.annotation.gtf.gz), optionally a custom catalogue (BED6/BED12).
Output: directory --out with the files:
  sites.tsv           all "sites" to be scored in step 2 (candidates, internal promoters, random controls)
  candidates_tx.tsv   candidates at the transcript level (with host gene information)
  host_genes.tsv      host genes (reference transcript, number of introns)
  host_tss.tsv        all annotated host start sites (for choosing the start site used in K562)
  sites.bed           the same sites as BED6 (for a genome browser)

Definitions (defaults, adjustable with parameters):
  host        = protein_coding gene
  host exons  = union of exons of host transcripts of type protein_coding / NMD / non_stop_decay /
                protein_coding_CDS_not_defined (NOT retained_intron, as these cover whole introns)
  candidate   = (a) an lncRNA transcript lying entirely within a host gene, on the same strand,
                not overlapping any host exon, or
                (b) a non-coding isoform of the host itself lying entirely within a reference intron
                (GENCODE often assigns such RNAs to the host gene); source: gencode_host_isoform
  context     = ref_intron if the candidate lies in an intron of the reference transcript
                (MANE Select > Ensembl canonical > longest), otherwise other_intronic
  internal promoter = start site of another host transcript lying in an intron of the reference
                transcript, at least --min-exon-dist from reference exons, where that transcript shares
                an exon with the reference (positive control for an "own start")
  random control = random position in reference introns of the same host, away from exons
                and from all known start sites ("pre-mRNA" control)

Coordinates in output files: 0-based (as in BED); `pos` = 0-based position of the start nucleotide.
"""
import argparse
import bisect
import gzip
import random
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')

HOST_EXON_TX_TYPES_DEFAULT = "protein_coding,nonsense_mediated_decay,non_stop_decay,protein_coding_CDS_not_defined"


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_gtf(path):
    tx = {}
    exons = defaultdict(list)
    n = 0
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "exon"):
                continue
            attrs, tags = {}, []
            for k, v in ATTR_RE.findall(f[8]):
                if k == "tag":
                    tags.append(v)
                else:
                    attrs.setdefault(k, v)
            tid = attrs.get("transcript_id")
            if f[2] == "transcript":
                tx[tid] = dict(
                    transcript_id=tid,
                    gene_id=attrs.get("gene_id"),
                    gene_name=attrs.get("gene_name", ""),
                    gene_type=attrs.get("gene_type", ""),
                    transcript_type=attrs.get("transcript_type", ""),
                    chrom=f[0], start=int(f[3]) - 1, end=int(f[4]), strand=f[6],
                    mane="MANE_Select" in tags, canonical="Ensembl_canonical" in tags,
                    readthrough="readthrough_transcript" in tags,
                )
                n += 1
            else:
                exons[tid].append((int(f[3]) - 1, int(f[4])))
    for tid in exons:
        exons[tid].sort()
    print(f"[GTF] transcripts: {n}", file=sys.stderr)
    return tx, exons


def parse_extra_bed(path):
    """Custom catalogue: BED6 or BED12. Returns a list of (transcript-like dict, exon list) pairs."""
    out = []
    with open_any(path) as fh:
        for i, line in enumerate(fh):
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                sys.exit(f"[ERROR] {path}: line {i+1} has fewer than 6 columns (stranded BED6/BED12 required)")
            chrom, s, e, name, strand = f[0], int(f[1]), int(f[2]), f[3], f[5]
            if len(f) >= 12:
                sizes = [int(x) for x in f[10].strip(",").split(",")]
                starts = [int(x) for x in f[11].strip(",").split(",")]
                ex = sorted((s + st, s + st + sz) for st, sz in zip(starts, sizes))
            else:
                ex = [(s, e)]
            out.append((dict(transcript_id=name, gene_id=name, gene_name=name, gene_type="extra",
                             transcript_type="extra", chrom=chrom, start=s, end=e, strand=strand,
                             mane=False, canonical=False), ex))
    print(f"[EXTRA] records: {len(out)}", file=sys.stderr)
    return out


def merge_intervals(ivs):
    ivs = sorted(ivs)
    out = []
    for s, e in ivs:
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


class IntervalSet:
    """Sorted, disjoint intervals [s, e) with fast overlap and distance queries."""

    def __init__(self, ivs):
        self.ivs = merge_intervals(ivs)
        self.starts = [s for s, _ in self.ivs]
        self.ends = [e for _, e in self.ivs]

    def overlaps(self, s, e):
        j = bisect.bisect_right(self.ends, s)
        return j < len(self.ivs) and self.starts[j] < e

    def dist(self, p):
        """Distance from position p to the nearest interval (0 if p lies inside one)."""
        if not self.ivs:
            return float("inf")
        j = bisect.bisect_right(self.starts, p)
        best = float("inf")
        if j > 0:
            s, e = self.ivs[j - 1]
            best = 0 if p < e else p - (e - 1)
        if j < len(self.ivs):
            best = min(best, self.starts[j] - p)
        return best


def tss_of(t):
    return t["start"] if t["strand"] == "+" else t["end"] - 1


def tes_of(t):
    return t["end"] - 1 if t["strand"] == "+" else t["start"]


def signed_dist(pos, ref, strand):
    """Distance along the direction of transcription (positive = downstream of the host start site)."""
    return pos - ref if strand == "+" else ref - pos


def exonic_len(ex):
    return sum(e - s for s, e in ex)


def pick_ref(tl, exons, host_exon_types):
    for key in ("mane", "canonical"):
        c = [t for t in tl if t[key]]
        if c:
            return c[0]
    c = [t for t in tl if t["transcript_type"] in host_exon_types] or tl
    return max(c, key=lambda t: exonic_len(exons[t["transcript_id"]]))


def build_hosts(tx, exons, host_gene_types, host_exon_types):
    by_gene = defaultdict(list)
    for t in tx.values():
        if t["gene_type"] in host_gene_types:
            by_gene[t["gene_id"]].append(t)
    hosts = {}
    for gid, tl in by_gene.items():
        if any(t.get("readthrough") for t in tl):
            continue  # readthrough genes (fusions of neighbouring genes) are not hosts
        typed = [t for t in tl if t["transcript_type"] in host_exon_types]
        if not typed:
            continue
        ref = pick_ref(tl, exons, host_exon_types)
        rex = exons[ref["transcript_id"]]
        introns = [(rex[i][1], rex[i + 1][0]) for i in range(len(rex) - 1) if rex[i + 1][0] > rex[i][1]]
        intron_starts = [s for s, _ in introns]
        strand = ref["strand"]

        def in_ref_intron(t):
            j = bisect.bisect_right(intron_starts, t["start"]) - 1
            return j >= 0 and introns[j][0] <= t["start"] and t["end"] <= introns[j][1]

        # host isoforms lying entirely within a reference intron (GENCODE often assigns sense intronic
        # RNAs to the host gene) -> these are candidates, not host exons
        intronic_iso = [t for t in tl if t["transcript_id"] != ref["transcript_id"]
                        and t["transcript_type"] != "protein_coding" and in_ref_intron(t)]
        iso_ids = {t["transcript_id"] for t in intronic_iso}
        typed_main = [t for t in typed if t["transcript_id"] not in iso_ids]
        if not typed_main:
            continue
        tss_up = min(tss_of(t) for t in typed_main) if strand == "+" else max(tss_of(t) for t in typed_main)
        hosts[gid] = dict(
            gene_id=gid, gene_name=ref["gene_name"], chrom=ref["chrom"], strand=strand,
            start=min(t["start"] for t in tl), end=max(t["end"] for t in tl),
            ref_tx=ref["transcript_id"], ref_tss=tss_of(ref), up_tss=tss_up,
            introns=introns,
            intron_starts=intron_starts,
            exon_set=IntervalSet([e for t in typed_main for e in exons[t["transcript_id"]]]),
            ref_exon_set=IntervalSet(rex),
            all_exon_set=IntervalSet([e for t in tl for e in exons[t["transcript_id"]]]),
            transcripts=typed_main,
            intronic_iso=intronic_iso,
        )
    print(f"[HOST] host genes: {len(hosts)}", file=sys.stderr)
    return hosts


def index_hosts(hosts):
    idx = defaultdict(list)
    for h in hosts.values():
        idx[(h["chrom"], h["strand"])].append(h)
    out = {}
    for k, hl in idx.items():
        hl.sort(key=lambda h: h["start"])
        maxlen = max(h["end"] - h["start"] for h in hl)
        out[k] = (hl, [h["start"] for h in hl], maxlen)
    return out


def containing_hosts(hidx, chrom, strand, s, e):
    if (chrom, strand) not in hidx:
        return []
    hl, starts, maxlen = hidx[(chrom, strand)]
    j0 = bisect.bisect_left(starts, s - maxlen)
    j = bisect.bisect_right(starts, s)
    return [h for h in hl[j0:j] if h["end"] >= e]


def ref_intron_of(h, s, e):
    j = bisect.bisect_right(h["intron_starts"], s) - 1
    if j >= 0:
        a, b = h["introns"][j]
        if a <= s and e <= b:
            n = len(h["introns"])
            num = j + 1 if h["strand"] == "+" else n - j
            return num, a, b
    return None


def site_record(site_id, site_type, h, pos, extra, exon_set=None):
    """exon_set: exons used for the distance (default: union of host exons;
    for internal promoters: reference exons, because the promoter's own first exon lies in an intron)."""
    intr = ref_intron_of(h, pos, pos + 1)
    exon_set = exon_set or h["exon_set"]
    rec = dict(
        site_id=site_id, site_type=site_type, chrom=h["chrom"], strand=h["strand"], pos=pos,
        host_gene_id=h["gene_id"], host_gene_name=h["gene_name"], host_ref_tx=h["ref_tx"],
        dist_ref_tss=signed_dist(pos, h["ref_tss"], h["strand"]),
        dist_up_tss=signed_dist(pos, h["up_tss"], h["strand"]),
        min_dist_host_exon=exon_set.dist(pos),
        intron_no=intr[0] if intr else "NA",
        intron_start=intr[1] if intr else "NA",
        intron_end=intr[2] if intr else "NA",
    )
    rec.update(extra)
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--extra-bed", help="custom iRNA catalogue (stranded BED6/BED12)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--host-gene-types", default="protein_coding")
    ap.add_argument("--host-exon-tx-types", default=HOST_EXON_TX_TYPES_DEFAULT)
    ap.add_argument("--cand-tx-types", default="lncRNA")
    ap.add_argument("--tss-cluster", type=int, default=50, help="distance for merging start sites of one gene (nt)")
    ap.add_argument("--min-exon-dist", type=int, default=1000, help="for internal promoters and controls (nt)")
    ap.add_argument("--ctrl-per-site", type=int, default=3)
    ap.add_argument("--ctrl-exclude", type=int, default=2000, help="exclusion zone around known start sites (nt)")
    ap.add_argument("--seed", type=int, default=20261004)
    ap.add_argument("--chrom", help="restrict to one chromosome (for testing, e.g. chr22)")
    a = ap.parse_args()

    import os
    os.makedirs(a.out, exist_ok=True)
    host_gene_types = set(a.host_gene_types.split(","))
    host_exon_types = set(a.host_exon_tx_types.split(","))
    cand_types = set(a.cand_tx_types.split(","))

    tx, exons = parse_gtf(a.gtf)
    if a.chrom:
        tx = {k: v for k, v in tx.items() if v["chrom"] == a.chrom}
    hosts = build_hosts(tx, exons, host_gene_types, host_exon_types)
    hidx = index_hosts(hosts)

    cand_items = [(t, exons[t["transcript_id"]]) for t in tx.values()
                  if t["transcript_type"] in cand_types and t["gene_type"] not in host_gene_types]
    if a.extra_bed:
        extra = parse_extra_bed(a.extra_bed)
        if a.chrom:
            extra = [x for x in extra if x[0]["chrom"] == a.chrom]
        cand_items += extra

    # --- candidates (transcript level) ---
    cand_rows = []
    n_multi = 0
    for t, ex in cand_items:
        hs = containing_hosts(hidx, t["chrom"], t["strand"], t["start"], t["end"])
        hs = [h for h in hs if not h["exon_set"].overlaps(t["start"], t["end"])]
        if not hs:
            continue
        if len(hs) > 1:
            n_multi += 1
        h = min(hs, key=lambda x: x["end"] - x["start"])  # innermost host
        intr = ref_intron_of(h, t["start"], t["end"])
        cand_rows.append(dict(
            cand_tx_id=t["transcript_id"], cand_gene_id=t["gene_id"], cand_gene_name=t["gene_name"],
            source="extra" if t["gene_type"] == "extra" else "gencode",
            chrom=t["chrom"], start=t["start"], end=t["end"], strand=t["strand"],
            tss=tss_of(t), tes=tes_of(t), n_exons=len(ex), exonic_len=exonic_len(ex),
            host_gene_id=h["gene_id"], host_gene_name=h["gene_name"], n_hosts=len(hs),
            context="ref_intron" if intr else "other_intronic",
            overlaps_any_host_exon_alltypes=h["all_exon_set"].overlaps(t["start"], t["end"]),
        ))
    n_sep = len(cand_rows)

    n_iso = 0
    for h in hosts.values():
        iso = sorted(h["intronic_iso"], key=lambda t: t["start"])
        loci = []
        for t in iso:
            if loci and t["start"] < loci[-1]["end"]:
                loci[-1]["end"] = max(loci[-1]["end"], t["end"])
                loci[-1]["tx"].append(t)
            else:
                loci.append(dict(end=t["end"], tx=[t]))
        for k, lc in enumerate(loci, 1):
            gid = f"{h['gene_id']}:iso{k}"
            for t in lc["tx"]:
                ex = exons[t["transcript_id"]]
                cand_rows.append(dict(
                    cand_tx_id=t["transcript_id"], cand_gene_id=gid, cand_gene_name=f"{h['gene_name']}-iso{k}",
                    source="gencode_host_isoform",
                    chrom=t["chrom"], start=t["start"], end=t["end"], strand=t["strand"],
                    tss=tss_of(t), tes=tes_of(t), n_exons=len(ex), exonic_len=exonic_len(ex),
                    host_gene_id=h["gene_id"], host_gene_name=h["gene_name"], n_hosts=1,
                    context="ref_intron", overlaps_any_host_exon_alltypes=True,
                ))
                n_iso += 1
    print(f"[CAND] candidate transcripts: {len(cand_rows)} (separate lncRNA/catalogue genes: {n_sep}, "
          f"host isoforms within an intron: {n_iso}; in >1 host: {n_multi})", file=sys.stderr)

    # --- candidate start units (merging start sites within a gene) ---
    by_gene = defaultdict(list)
    for r in cand_rows:
        by_gene[r["cand_gene_id"]].append(r)
    sites = []
    known_tss = defaultdict(list)  # (chrom,strand) -> positions of known start sites (excluded zones for controls)
    for gid, rows in by_gene.items():
        rows.sort(key=lambda r: r["tss"])
        clusters = []
        for r in rows:
            if clusters and abs(r["tss"] - clusters[-1][-1]["tss"]) <= a.tss_cluster:
                clusters[-1].append(r)
            else:
                clusters.append([r])
        for cl in clusters:
            # unit position = most frequent start site (on ties: the most upstream one)
            counts = defaultdict(int)
            for r in cl:
                counts[r["tss"]] += 1
            best = max(counts.items(), key=lambda kv: (kv[1], -kv[0] if cl[0]["strand"] == "+" else kv[0]))[0]
            h = hosts[cl[0]["host_gene_id"]]
            sid = f"{gid}:{best}"
            sites.append(site_record(sid, "candidate", h, best, dict(
                cand_gene_id=gid, cand_gene_name=cl[0]["cand_gene_name"], source=cl[0]["source"],
                cand_tx_ids=",".join(sorted({r["cand_tx_id"] for r in cl})),
                tes_list=",".join(str(x) for x in sorted({r["tes"] for r in cl})),
                context=",".join(sorted({r["context"] for r in cl})),
            )))
            known_tss[(h["chrom"], h["strand"])].append(best)

    # --- internal promoters of the host (positive control) ---
    n_int = 0
    for h in hosts.values():
        seen = []
        for t in h["transcripts"]:
            if t["transcript_id"] == h["ref_tx"]:
                continue
            p = tss_of(t)
            if ref_intron_of(h, p, p + 1) is None:
                continue
            if h["ref_exon_set"].dist(p) < a.min_exon_dist:
                continue
            if not any(h["ref_exon_set"].overlaps(es, ee) for es, ee in exons[t["transcript_id"]]):
                continue  # shares no exon with the host reference -> not an alternative promoter
            if any(abs(p - q) <= a.tss_cluster for q in seen):
                continue
            seen.append(p)
            sites.append(site_record(f"{h['gene_id']}:int:{p}", "internal_promoter", h, p, dict(
                cand_gene_id="NA", cand_gene_name="NA", source="gencode",
                cand_tx_ids=t["transcript_id"], tes_list="NA", context="ref_intron"), exon_set=h["ref_exon_set"]))
            known_tss[(h["chrom"], h["strand"])].append(p)
            n_int += 1
    print(f"[CTRL+] internal promoters: {n_int}", file=sys.stderr)

    # --- random controls in introns of the same host ---
    for k in known_tss:
        known_tss[k].sort()
    rng = random.Random(a.seed)
    allowed_cache = {}

    def allowed_ivs(h):
        if h["gene_id"] in allowed_cache:
            return allowed_cache[h["gene_id"]]
        kt = known_tss[(h["chrom"], h["strand"])]
        ivs = []
        for s, e in h["introns"]:
            s2, e2 = s + a.min_exon_dist, e - a.min_exon_dist
            if e2 <= s2:
                continue
            # cut out zones around known start sites
            cut = [(max(s2, p - a.ctrl_exclude), min(e2, p + a.ctrl_exclude + 1))
                   for p in kt[bisect.bisect_left(kt, s2 - a.ctrl_exclude):bisect.bisect_right(kt, e2 + a.ctrl_exclude)]]
            cur = s2
            for cs, ce in merge_intervals([c for c in cut if c[1] > c[0]]):
                if cs > cur:
                    ivs.append((cur, cs))
                cur = max(cur, ce)
            if cur < e2:
                ivs.append((cur, e2))
        # the distance from exons of other host transcripts is checked when sampling
        allowed_cache[h["gene_id"]] = ivs
        return ivs

    n_ctrl = 0
    cand_sites = [s for s in sites if s["site_type"] == "candidate"]
    for s in cand_sites:
        h = hosts[s["host_gene_id"]]
        ivs = allowed_ivs(h)
        tot = sum(e - b for b, e in ivs)
        if tot <= 0:
            continue
        for i in range(a.ctrl_per_site):
            for _try in range(20):
                x = rng.randrange(tot)
                for b, e in ivs:
                    if x < e - b:
                        p = b + x
                        break
                    x -= e - b
                if h["exon_set"].dist(p) >= a.min_exon_dist:
                    break
            else:
                continue
            sites.append(site_record(f"{s['site_id']}:ctrl{i+1}", "random_control", h, p, dict(
                cand_gene_id=s["cand_gene_id"], cand_gene_name="NA", source="random",
                cand_tx_ids="NA", tes_list="NA", context="ref_intron")))
            n_ctrl += 1
    print(f"[CTRL-] random controls: {n_ctrl}", file=sys.stderr)

    # --- output ---
    cols = ["site_id", "site_type", "chrom", "strand", "pos", "host_gene_id", "host_gene_name", "host_ref_tx",
            "dist_ref_tss", "dist_up_tss", "min_dist_host_exon", "intron_no", "intron_start", "intron_end",
            "cand_gene_id", "cand_gene_name", "source", "cand_tx_ids", "tes_list", "context"]
    with open(f"{a.out}/sites.tsv", "w") as fo, open(f"{a.out}/sites.bed", "w") as fb:
        fo.write("\t".join(cols) + "\n")
        for r in sorted(sites, key=lambda r: (r["chrom"], r["pos"])):
            fo.write("\t".join(str(r[c]) for c in cols) + "\n")
            fb.write(f"{r['chrom']}\t{r['pos']}\t{r['pos']+1}\t{r['site_id']}|{r['site_type']}\t0\t{r['strand']}\n")

    ccols = list(cand_rows[0].keys()) if cand_rows else []
    with open(f"{a.out}/candidates_tx.tsv", "w") as fo:
        fo.write("\t".join(ccols) + "\n")
        for r in cand_rows:
            fo.write("\t".join(str(r[c]) for c in ccols) + "\n")

    with open(f"{a.out}/host_genes.tsv", "w") as fo, open(f"{a.out}/host_tss.tsv", "w") as ft:
        fo.write("host_gene_id\thost_gene_name\tchrom\tstrand\tstart\tend\tref_tx\tref_tss\tup_tss\tn_ref_introns\n")
        ft.write("host_gene_id\ttranscript_id\tchrom\tstrand\tpos\tis_ref\n")
        for h in hosts.values():
            fo.write(f"{h['gene_id']}\t{h['gene_name']}\t{h['chrom']}\t{h['strand']}\t{h['start']}\t{h['end']}\t"
                     f"{h['ref_tx']}\t{h['ref_tss']}\t{h['up_tss']}\t{len(h['introns'])}\n")
            for t in h["transcripts"]:
                ft.write(f"{h['gene_id']}\t{t['transcript_id']}\t{h['chrom']}\t{h['strand']}\t{tss_of(t)}\t"
                         f"{int(t['transcript_id'] == h['ref_tx'])}\n")

    n_cand_genes = len({s["cand_gene_id"] for s in cand_sites})
    print(f"[OUT] candidate genes: {n_cand_genes}; start units: {len(cand_sites)}; written to {a.out}/",
          file=sys.stderr)


if __name__ == "__main__":
    main()
