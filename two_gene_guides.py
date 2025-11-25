#!/usr/bin/env python3
import argparse, os, sys, csv, re, json
from typing import Dict, List, Tuple, Any, Optional
import gffutils
from pymongo import MongoClient
from pprint import pprint

# -----------------------------
# Utility: CSV writers & logging
# -----------------------------
def write_csv(rows: List[dict], path: str, field_order: List[str] = None, max_preview: int = 5):
    if not rows:
        if field_order:
            with open(path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=field_order)
                w.writeheader()
        print(f"[write_csv] No rows; created empty file with header: {path}")
        return
    if field_order is None:
        keys = set()
        for r in rows: keys.update(r.keys())
        field_order = list(sorted(keys))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=field_order)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in field_order})
    print(f"[write_csv] Wrote {len(rows)} rows -> {path}")
    print("[write_csv] Preview:")
    for r in rows[:max_preview]:
        print("  ", {k: r.get(k, "") for k in field_order})

def log_header(title: str):
    print("\n" + "="*80)
    print(title)
    print("="*80)

# -----------------------------
# Chromosome normalization
# -----------------------------
def normalize_chrom_name(chrom: str, prefix: str = "chr") -> str:
    """
    Convert GFF chrom names (e.g., 'Chr4_RagTag', 'chr04', '4') to Mongo style like 'chr4'.
    Strips suffixes after underscore/dash/dot/space, drops leading zeros, lowercases.
    """
    if not chrom:
        return chrom
    base = re.split(r'[_\-\.\s]', str(chrom))[0]  # 'Chr4_RagTag' -> 'Chr4'
    m_num = re.match(r'^(?:chr|Chr)?0*([0-9]+)$', base)
    if m_num:
        return f"{prefix}{int(m_num.group(1))}"
    m_num_prefix = re.match(r'^(?:chr|Chr)?0*([0-9]+)', base)
    if m_num_prefix:
        return f"{prefix}{int(m_num_prefix.group(1))}"
    m_letters = re.match(r'^(?:chr|Chr)?([A-Za-z]+)$', base)
    if m_letters:
        return f"{prefix}{m_letters.group(1).lower()}"
    return base.lower()

def normalize_intervals_chrom(intervals, prefix: str = "chr", verbose: bool = True):
    out = []
    for iv in intervals:
        old = iv["chrom"]
        new = normalize_chrom_name(old, prefix=prefix)
        if verbose and old != new:
            print(f"[ChromNormalize] {old} -> {new}")
        out.append({**iv, "chrom": new})
    return out

# -----------------------------
# GFF database helpers
# -----------------------------
def open_or_build_gffdb(gff_db_path: str, gff_fallback_path: str = None) -> gffutils.FeatureDB:
    if os.path.exists(gff_db_path):
        print(f"[GFF] Opening DB: {gff_db_path}")
        return gffutils.FeatureDB(gff_db_path, keep_order=True)
    if not gff_fallback_path:
        raise FileNotFoundError(f"GFF DB not found: {gff_db_path} (no fallback GFF3 provided to build from)")
    print(f"[GFF] Building DB from {gff_fallback_path} -> {gff_db_path}")
    gffutils.create_db(
        gff_fallback_path,
        dbfn=gff_db_path,
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        force_gff=True
    )
    return gffutils.FeatureDB(gff_db_path, keep_order=True)

def resolve_target(db: gffutils.FeatureDB, id_str: str):
    """
    Resolve an input ID to:
      - parent gene feature,
      - list of transcript IDs to use (restrict to one if a transcript was provided),
      - a label describing the resolution.
    """
    # 1) direct
    try:
        feat = db[id_str]
        if feat.featuretype in ("mRNA", "transcript"):
            # parent gene
            parents = list(db.parents(feat, featuretype="gene"))
            if not parents:
                raise gffutils.FeatureNotFoundError(f"No parent gene found for transcript {id_str}")
            tx_ids = [feat.id]
            return parents[0], tx_ids, f"transcript:{feat.id}"
        if feat.featuretype == "gene":
            # all transcripts under the gene
            tx_ids = [tx.id for tx in db.children(feat, featuretype=("mRNA","transcript"))]
            return feat, tx_ids, f"gene:{feat.id}"
    except Exception:
        pass

    # 2) fallback: scan genes by attributes
    for g in db.features_of_type("gene"):
        a = g.attributes
        if id_str in (a.get("ID",[None])[0], a.get("Name",[None])[0],
                      a.get("gene_id",[None])[0], a.get("gene_name",[None])[0]):
            tx_ids = [tx.id for tx in db.children(g, featuretype=("mRNA","transcript"))]
            return g, tx_ids, f"gene_attr:{id_str}"

    # 3) heuristic: treat as transcript with suffix or gene base
    if "." in id_str:
        # try as exact transcript ID
        try:
            t = db[id_str]
            if t.featuretype in ("mRNA","transcript"):
                parents = list(db.parents(t, featuretype="gene"))
                if not parents:
                    raise gffutils.FeatureNotFoundError
                return parents[0], [t.id], f"transcript_heuristic:{t.id}"
        except Exception:
            # try base as gene
            base = id_str.split(".")[0]
            try:
                g = db[base]
                if g.featuretype == "gene":
                    tx_ids = [tx.id for tx in db.children(g, featuretype=("mRNA","transcript"))]
                    return g, tx_ids, f"gene_from_base:{base}"
            except Exception:
                pass

    raise gffutils.FeatureNotFoundError(f"Could not resolve '{id_str}' to gene/transcript")

def cds_intervals_per_transcript(db: gffutils.FeatureDB, gene_feat: gffutils.Feature,
                                 restrict_tx_ids: Optional[List[str]] = None):
    """
    Return dict mapping transcript_id -> list of CDS tuples (chrom, start, end, strand).
    If restrict_tx_ids is provided, only those transcripts are used.
    Also includes CDS directly under the gene grouped by Parent attribute.
    """
    from collections import defaultdict
    cds_by_tx = defaultdict(list)

    # transcripts under gene
    for tx in db.children(gene_feat, featuretype=("mRNA","transcript"), order_by="start"):
        if restrict_tx_ids and tx.id not in restrict_tx_ids:
            continue
        for cds in db.children(tx, featuretype="CDS", order_by="start"):
            cds_by_tx[tx.id].append((cds.chrom, cds.start, cds.end, cds.strand))

    # CDS directly under gene (group by Parent)
    for cds in db.children(gene_feat, featuretype="CDS", order_by="start"):
        parents = cds.attributes.get("Parent", [])
        for p in parents or [gene_feat.id]:
            if restrict_tx_ids and p not in (restrict_tx_ids + [gene_feat.id]):
                continue
            cds_by_tx[p].append((cds.chrom, cds.start, cds.end, cds.strand))

    return dict(cds_by_tx)

def merge_intervals_union(cds_map: Dict[str, List[Tuple[str,int,int,str]]]):
    all_parts, strands = [], set()
    for parts in cds_map.values():
        all_parts.extend(parts)
        for _,_,_,s in parts:
            strands.add(s)
    if not all_parts:
        return []

    all_parts.sort(key=lambda x: (x[0], x[1], x[2]))
    merged = []
    cur_chrom, cur_start, cur_end = all_parts[0][0], all_parts[0][1], all_parts[0][2]
    for chrom, s, e, _ in all_parts[1:]:
        if chrom != cur_chrom or s > cur_end + 1:
            merged.append({"chrom": cur_chrom, "start": cur_start, "end": cur_end})
            cur_chrom, cur_start, cur_end = chrom, s, e
        else:
            if e > cur_end:
                cur_end = e
    merged.append({"chrom": cur_chrom, "start": cur_start, "end": cur_end})
    strand = list(strands)[0] if len(strands) == 1 else "."
    for m in merged:
        m["strand"] = strand
    return merged

def write_bed(intervals: List[Dict[str, Any]], bed_path: str, name_prefix: str):
    with open(bed_path, "w") as out:
        for i, iv in enumerate(intervals, 1):
            out.write(f"{iv['chrom']}\t{iv['start']-1}\t{iv['end']}\t{name_prefix}_{i}\t0\t{iv.get('strand','+')}\n")

# -----------------------------
# Mongo helpers
# -----------------------------
def pam_filter_clause(pam: str) -> dict:
    if pam.upper() == "NGG":
        # exact triplets; uses the PAM_sequence equality in the index
        return {"$or": [
            {"pam_type": "NGG"},
            {"PAM_sequence": {"$in": ["TGG", "AGG", "GGG", "CGG"]}},
        ]}
    else:
        return {"$or": [
            {"pam_type": pam},
            {"PAM_sequence": pam}
        ]}

def intervals_to_or_clauses(intervals: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    ors = []
    for iv in intervals:
        ors.append({
            "$and": [
                {"chromosome": iv["chrom"]},
                {"$or": [
                    {"cutsite": {"$gte": iv["start"], "$lte": iv["end"]}},
                    {"$and": [  # spacer overlap fallback
                        {"start_pos": {"$lte": iv["end"]}},
                        {"end_pos":   {"$gte": iv["start"]}}
                    ]}
                ]}
            ]
        })
    return ors

def fetch_guides_in_intervals(col, intervals, pam, score_min=None, fields=None, require_full_span=False):
    if fields is None:
        fields = {"sequence":1,"PAM_sequence":1,"pam_type":1,"chromosome":1,"cutsite":1,"strand":1,
                  "on_site_score":1,"start_pos":1,"end_pos":1}
    or_part = intervals_fullspan_or(intervals) if require_full_span else intervals_to_or_clauses(intervals)
    filt = {"$and": [ {"$or": or_part}, pam_filter_clause(pam) ]}
    if score_min is not None:
        filt["on_site_score"] = {"$gte": score_min}
    # optional: hint so Mongo uses your index
    cur = col.find(filt, fields).hint([("chromosome",1),("cutsite",1),("PAM_sequence",1)])
    return list(cur)

def fetch_occurrence_count(counts_col, seq: str, pam: str) -> int:
    doc = counts_col.find_one({"seq": seq, "pam": pam}, {"occurrence_count": 1})
    if doc and "occurrence_count" in doc:
        return int(doc["occurrence_count"])
    doc = counts_col.find_one({"_id": {"seq": seq, "pam": pam}}, {"occurrence_count": 1})
    if doc and "occurrence_count" in doc:
        return int(doc["occurrence_count"])
    return 0

def fetch_all_loci_for_guide(col, seq: str, pam: str, fields=None):
    if fields is None:
        fields = {"sequence":1, "PAM_sequence":1, "chromosome":1, "cutsite":1, "strand":1, "on_site_score":1, "start_pos":1, "end_pos":1}
    cur = col.find({"sequence": seq, "PAM_sequence": pam}, fields)
    seen, loci = set(), []
    for d in cur:
        key = (d.get("chromosome"), d.get("cutsite"), d.get("strand"))
        if key not in seen:
            seen.add(key)
            loci.append(d)
    return loci

def overlaps_iv(chrom: str, pos: int, ivs: List[Dict[str, Any]]):
    chrom = normalize_chrom_name(chrom)
    for iv in ivs:
        if chrom == iv["chrom"] and iv["start"] <= pos <= iv["end"]:
            return True
    return False

# Require full span of intervals (i.e., guide fully inside CDS)

def intervals_fullspan_or(intervals):
    # start_pos >= CDS_start AND end_pos <= CDS_end
    return [{"$and": [
              {"chromosome": iv["chrom"]},
              {"start_pos": {"$gte": iv["start"]}},
              {"end_pos":   {"$lte": iv["end"]}}
            ]} for iv in intervals]


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Verbose pipeline to find guides that hit both genes' (or specified transcripts') CDS with occurrence_count==2."
    )
    ap.add_argument("gene_name")
    ap.add_argument("gene_id_1", help="Gene ID OR transcript ID (e.g., g15664 OR g15664.t2)")
    ap.add_argument("gene_id_2", help="Gene ID OR transcript ID (e.g., g36136 OR g36136.t1)")
    ap.add_argument("--gffdb", default="braker.gff.db")
    ap.add_argument("--gff", default=None)
    ap.add_argument("--mongo-uri", dest="mongo_uri", default="mongodb://127.0.0.1:27017")
    ap.add_argument("--db", default="CROPSR")
    ap.add_argument("--guides-col", dest="guides_col", default="guides")
    ap.add_argument("--counts-col", dest="counts_col", default="guide_counts")
    ap.add_argument("--pam", default="NGG")
    ap.add_argument("--score-min", dest="score_min", type=float, default=0.10)
    ap.add_argument("--write-bed", dest="write_bed", action="store_true")
    ap.add_argument("--out", default=None)

    args = ap.parse_args()
    out_csv = args.out or f"{args.gene_name}_validated_guides.csv"

    # 1) GFF: open DB & resolve targets
    log_header("STEP 1: Open GFF DB and resolve input IDs (gene or transcript)")
    db = open_or_build_gffdb(args.gffdb, args.gff)
    gene1, txs1, label1 = resolve_target(db, args.gene_id_1)
    gene2, txs2, label2 = resolve_target(db, args.gene_id_2)
    print(f"[GFF] Target1: {label1} -> gene {gene1.id} {gene1.chrom}:{gene1.start}-{gene1.end} {gene1.strand}")
    print(f"[GFF] Target2: {label2} -> gene {gene2.id} {gene2.chrom}:{gene2.start}-{gene2.end} {gene2.strand}")
    print(f"[GFF] Using transcripts:")
    print(f"  Target1 transcripts: {txs1}")
    print(f"  Target2 transcripts: {txs2}")

    # 2) Extract CDS per (restricted) transcript & union
    log_header("STEP 2: Extract CDS intervals (restricted to requested transcripts if provided)")
    cds1_map = cds_intervals_per_transcript(db, gene1, restrict_tx_ids=txs1)
    cds2_map = cds_intervals_per_transcript(db, gene2, restrict_tx_ids=txs2)
    print(f"[GFF] Target1 transcripts with CDS (after restriction): {len(cds1_map)}")
    print(f"[GFF] Target2 transcripts with CDS (after restriction): {len(cds2_map)}")
    for tx, parts in list(cds1_map.items())[:2]:
        print("[GFF] Target1 example TX:", tx, "CDS parts:", len(parts), "first:", parts[:2])
    for tx, parts in list(cds2_map.items())[:2]:
        print("[GFF] Target2 example TX:", tx, "CDS parts:", len(parts), "first:", parts[:2])

    cds1_union = merge_intervals_union(cds1_map)
    cds2_union = merge_intervals_union(cds2_map)
    cds1_union = normalize_intervals_chrom(cds1_union, prefix="chr", verbose=True)
    cds2_union = normalize_intervals_chrom(cds2_union, prefix="chr", verbose=True)
    print("[GFF] Target1 CDS union intervals:", len(cds1_union))
    print("[GFF] Target2 CDS union intervals:", len(cds2_union))
    if args.write_bed:
        write_bed(cds1_union, f"{args.gene_name}_{args.gene_id_1}_CDS_union.bed", f"{args.gene_id_1}_CDS")
        write_bed(cds2_union, f"{args.gene_name}_{args.gene_id_2}_CDS_union.bed", f"{args.gene_id_2}_CDS")
        print("[BED] Wrote CDS union BEDs.")

    # 3) Mongo: connect
    log_header("STEP 3: Connect to MongoDB")
    client = MongoClient(args.mongo_uri)
    G = client[args.db][args.guides_col]
    GC = client[args.db][args.counts_col]
    print("[Mongo] Connected. Collections:", args.guides_col, args.counts_col)

    # 4) Query guides within CDS unions
    log_header("STEP 4: Query guides within CDS union for each target")
    guides1 = fetch_guides_in_intervals(G, cds1_union, args.pam, args.score_min)
    print(f"[Mongo] Target1 guides in CDS: {len(guides1)}")
    if guides1[:3]: pprint(guides1[:3])
    guides2 = fetch_guides_in_intervals(G, cds2_union, args.pam, args.score_min)
    print(f"[Mongo] Target2 guides in CDS: {len(guides2)}")
    if guides2[:3]: pprint(guides2[:3])

    write_csv(guides1, f"{args.gene_name}_target1_guides_raw.csv")
    write_csv(guides2, f"{args.gene_name}_target2_guides_raw.csv")

    # 5) Intersect by (sequence, PAM_sequence)
    log_header("STEP 5: Intersect guides by (sequence, PAM)")
    def key_of(d): return (d.get("sequence"), d.get("PAM_sequence"))
    set1 = { key_of(d) for d in guides1 if d.get("sequence") and d.get("PAM_sequence") }
    set2 = { key_of(d) for d in guides2 if d.get("sequence") and d.get("PAM_sequence") }
    common_keys = sorted(list(set1 & set2))
    print(f"[Intersect] common (sequence,PAM) count = {len(common_keys)}")
    if common_keys[:10]: print("[Intersect] sample keys:", common_keys[:10])
    write_csv([{"sequence": s, "PAM_sequence": p} for (s,p) in common_keys],
              f"{args.gene_name}_common_keys.csv", ["sequence","PAM_sequence"])

    # 6) occurrence_count==2 filter
    log_header("STEP 6: Check occurrence_count==2 in guide_counts for common keys")
    occ2_keys, occ_rows = [], []
    for (seq, pam) in common_keys:
        occ = fetch_occurrence_count(GC, seq, pam)
        occ_rows.append({"sequence": seq, "PAM_sequence": pam, "occurrence_count": occ})
        if occ == 2:
            occ2_keys.append((seq, pam))
    print(f"[Counts] common keys with occ==2: {len(occ2_keys)} / {len(common_keys)}")
    write_csv(occ_rows, f"{args.gene_name}_common_keys_with_counts.csv",
              ["sequence","PAM_sequence","occurrence_count"])

    # 7) Validate loci belong one to each target's CDS
    log_header("STEP 7: Validate two loci map one to each target's CDS")
    accepted_rows, debug_rows = [], []
    for (seq, pam) in occ2_keys:
        loci = fetch_all_loci_for_guide(G, seq, pam)
        for d in loci:
            debug_rows.append({
                "sequence": seq, "PAM_sequence": pam,
                "chromosome": d.get("chromosome"),
                "cutsite": d.get("cutsite"),
                "strand": d.get("strand"),
                "on_site_score": d.get("on_site_score")
            })
        if len(loci) != 2:
            continue
        inside1 = [d for d in loci if overlaps_iv(d.get("chromosome"), int(d.get("cutsite")), cds1_union)]
        inside2 = [d for d in loci if overlaps_iv(d.get("chromosome"), int(d.get("cutsite")), cds2_union)]
        if not inside1 or not inside2:
            continue
        d1, d2 = inside1[0], inside2[0]
        row = {
            "gene_name": args.gene_name,
            "gene_id_1": args.gene_id_1,
            "gene_id_2": args.gene_id_2,
            "sequence": seq,
            "PAM_sequence": pam,
            "occurrence_count": 2,
            "g1_chr": d1.get("chromosome"),
            "g1_cutsite": d1.get("cutsite"),
            "g1_strand": d1.get("strand"),
            "g1_on_site_score": d1.get("on_site_score"),
            "g2_chr": d2.get("chromosome"),
            "g2_cutsite": d2.get("cutsite"),
            "g2_strand": d2.get("strand"),
            "g2_on_site_score": d2.get("on_site_score"),
        }
        accepted_rows.append(row)

    print(f"[Validate] accepted guides: {len(accepted_rows)}")
    write_csv(debug_rows, f"{args.gene_name}_occ2_loci_debug.csv",
              ["sequence","PAM_sequence","chromosome","cutsite","strand","on_site_score"])

    # 8) Final CSV
    log_header("STEP 8: Write final validated guides CSV")
    if accepted_rows:
        field_order = [
            "gene_name","gene_id_1","gene_id_2",
            "sequence","PAM_sequence","occurrence_count",
            "g1_chr","g1_cutsite","g1_strand","g1_on_site_score",
            "g2_chr","g2_cutsite","g2_strand","g2_on_site_score"
        ]
        write_csv(accepted_rows, (args.out or f"{args.gene_name}_validated_guides.csv"), field_order)
    else:
        print("[Final] No validated guides found. Inspect intermediate CSVs:")
        print(f"  - {args.gene_name}_target1_guides_raw.csv")
        print(f"  - {args.gene_name}_target2_guides_raw.csv")
        print(f"  - {args.gene_name}_common_keys.csv")
        print(f"  - {args.gene_name}_common_keys_with_counts.csv")
        print(f"  - {args.gene_name}_occ2_loci_debug.csv")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR:", e)
        sys.exit(1)


