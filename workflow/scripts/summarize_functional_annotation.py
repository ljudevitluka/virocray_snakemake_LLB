#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from datetime import datetime

import pandas as pd

IPS_COLUMNS = [
    "orf_id", "md5", "orf_length", "analysis", "signature_acc", "signature_desc",
    "start", "stop", "score", "status", "date", "interpro_acc", "interpro_desc",
]
PROFILE_ANALYSES = {"ProSiteProfiles", "ProSitePatterns", "HAMAP"}


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarise InterProScan and NeoRdRp2 results per representative contig, including RdRp calls"
    )
    p.add_argument("--reps", required=True, help="Representative contigs FASTA (headers <contig_id>_cluster_<n>)")
    p.add_argument("--orfs", required=True, help="getorf protein FASTA (headers <rep header>_<n>)")
    p.add_argument("--interproscan", required=True, help="InterProScan TSV output")
    p.add_argument("--neordrp", required=True, nargs="+", help="hmmsearch --domtblout files against NeoRdRp2")
    p.add_argument("--rdrp-interpro-ids", required=True, help="TSV with column interpro_id listing RdRp InterPro entries")
    p.add_argument("--rdrp-signature-regex", required=True, help="Pattern for RdRp signatures without an InterPro entry")
    p.add_argument("--exclude-signatures", nargs="*", default=[], help="Signature accessions never counted as RdRp")
    p.add_argument("--neordrp-evalue", type=float, required=True, help="Full-sequence E-value cutoff for NeoRdRp2 hits")
    p.add_argument("--out-table", required=True)
    p.add_argument("--out-rdrp-faa", required=True, help="FASTA of the best RdRp ORF per RdRp-positive contig")
    return p.parse_args()


def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def rep_to_contig(rep_id):
    return rep_id.rsplit("_cluster_", 1)[0]


def orf_to_contig(orf_id):
    # getorf appends _<n> to the representative header
    return rep_to_contig(orf_id.rsplit("_", 1)[0])


def read_interproscan(path):
    try:
        df = pd.read_csv(path, sep="\t", header=None, dtype=str, usecols=range(13), names=IPS_COLUMNS)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=IPS_COLUMNS, dtype=str)
    df["contig_id"] = df["orf_id"].map(orf_to_contig)
    # Column 9 is an E-value for HMM-based databases but a match score for profile-based ones
    score = pd.to_numeric(df["score"], errors="coerce")
    is_profile = df["analysis"].isin(PROFILE_ANALYSES)
    df["evalue"] = score.where(~is_profile)
    df["profile_score"] = score.where(is_profile)
    return df


def read_domtbl(paths, max_evalue):
    rows = []
    for path in paths:
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split()
                rows.append({
                    "orf_id": c[0], "hmm": c[3], "hmm_len": int(c[5]),
                    "evalue": float(c[6]), "bitscore": float(c[7]),
                    "hmm_from": int(c[15]), "hmm_to": int(c[16]),
                })
    df = pd.DataFrame(rows, columns=["orf_id", "hmm", "hmm_len", "evalue", "bitscore", "hmm_from", "hmm_to"])
    df = df[df["evalue"] <= max_evalue]
    df["contig_id"] = df["orf_id"].map(orf_to_contig)
    return df


def hmm_coverage(group):
    # Fraction of the HMM covered by the union of its domain alignments on this ORF
    covered, end = 0, 0
    for s, e in sorted(zip(group["hmm_from"], group["hmm_to"])):
        if e > end:
            covered += e - max(s, end + 1) + 1
            end = e
    return covered / group["hmm_len"].iloc[0]


def summarize_neordrp(domtbl):
    cols = ["contig_id", "neordrp_best_orf", "neordrp_best_hmm", "neordrp_evalue",
            "neordrp_bitscore", "neordrp_hmm_coverage", "neordrp_n_hmms"]
    if domtbl.empty:
        return pd.DataFrame(columns=cols)

    pairs = (domtbl.groupby(["contig_id", "orf_id", "hmm"])
             .agg(evalue=("evalue", "first"), bitscore=("bitscore", "first"))
             .reset_index())
    cov = domtbl.groupby(["orf_id", "hmm"]).apply(hmm_coverage).rename("hmm_coverage").reset_index()
    pairs = pairs.merge(cov, on=["orf_id", "hmm"])

    best = pairs.sort_values(["evalue", "bitscore"], ascending=[True, False]).groupby("contig_id").first()
    n_hmms = pairs.groupby("contig_id")["hmm"].nunique()
    return pd.DataFrame({
        "contig_id": best.index,
        "neordrp_best_orf": best["orf_id"].values,
        "neordrp_best_hmm": best["hmm"].values,
        "neordrp_evalue": best["evalue"].values,
        "neordrp_bitscore": best["bitscore"].values,
        "neordrp_hmm_coverage": best["hmm_coverage"].round(3).values,
        "neordrp_n_hmms": n_hmms.loc[best.index].values,
    })


def flag_rdrp_matches(ips, rdrp_ids, regex, exclude):
    unintegrated = ips["interpro_acc"].isna() | (ips["interpro_acc"] == "-")
    by_entry = ips["interpro_acc"].isin(rdrp_ids)
    by_name = unintegrated & ips["signature_desc"].fillna("").str.contains(regex, flags=re.IGNORECASE, regex=True)
    return (by_entry | by_name) & ~ips["signature_acc"].isin(exclude)


def summarize_ips_rdrp(rdrp):
    cols = ["contig_id", "ips_rdrp_best_orf", "ips_rdrp_best_signature", "ips_rdrp_evalue",
            "ips_rdrp_profile_score", "ips_rdrp_interpro", "ips_rdrp_n_databases", "ips_rdrp_signatures"]
    if rdrp.empty:
        return pd.DataFrame(columns=cols)

    rdrp = rdrp.copy()
    rdrp["label"] = rdrp["analysis"] + ":" + rdrp["signature_acc"] + " (" + rdrp["signature_desc"].fillna("-") + ")"
    out = []
    for contig, g in rdrp.groupby("contig_id"):
        # Best = lowest E-value; profile-based matches (no E-value) sort last
        best = g.sort_values("evalue", na_position="last").iloc[0]
        ipr = sorted({a for a in g["interpro_acc"] if isinstance(a, str) and a != "-"})
        out.append({
            "contig_id": contig,
            "ips_rdrp_best_orf": best["orf_id"],
            "ips_rdrp_best_signature": best["label"],
            "ips_rdrp_evalue": best["evalue"],
            "ips_rdrp_profile_score": g["profile_score"].max(),
            "ips_rdrp_interpro": ";".join(ipr) if ipr else "-",
            "ips_rdrp_n_databases": g["analysis"].nunique(),
            "ips_rdrp_signatures": ";".join(sorted(set(g["label"]))),
        })
    return pd.DataFrame(out, columns=cols)


def summarize_general(ips):
    rows = []
    for contig, g in ips.groupby("contig_id"):
        ipr = (g[g["interpro_acc"].notna() & (g["interpro_acc"] != "-")]
               .drop_duplicates("interpro_acc").sort_values("interpro_acc"))
        pfam = g[g["analysis"] == "Pfam"].drop_duplicates("signature_acc").sort_values("signature_acc")
        rows.append({
            "contig_id": contig,
            "n_orfs_annotated": g["orf_id"].nunique(),
            "interpro_annotations": ";".join(ipr["interpro_acc"] + ":" + ipr["interpro_desc"].fillna("-")) or "-",
            "pfam_annotations": ";".join(pfam["signature_acc"] + ":" + pfam["signature_desc"].fillna("-")) or "-",
        })
    return pd.DataFrame(rows, columns=["contig_id", "n_orfs_annotated", "interpro_annotations", "pfam_annotations"])


def main():
    args = parse_args()

    reps = read_fasta(args.reps)
    orfs = read_fasta(args.orfs)
    table = pd.DataFrame({"contig_id": [rep_to_contig(r) for r in reps]})
    n_orfs = pd.Series([orf_to_contig(o) for o in orfs]).value_counts()
    table["n_orfs"] = table["contig_id"].map(n_orfs).fillna(0).astype(int)
    log(f"Representatives: {len(table)}, ORFs: {len(orfs)}")

    ips = read_interproscan(args.interproscan)
    rdrp_ids = set(pd.read_csv(args.rdrp_interpro_ids, sep="\t")["interpro_id"])
    ips_rdrp = ips[flag_rdrp_matches(ips, rdrp_ids, args.rdrp_signature_regex, set(args.exclude_signatures))]
    neordrp = read_domtbl(args.neordrp, args.neordrp_evalue)

    table = (table
             .merge(summarize_general(ips), on="contig_id", how="left")
             .merge(summarize_neordrp(neordrp), on="contig_id", how="left")
             .merge(summarize_ips_rdrp(ips_rdrp), on="contig_id", how="left"))
    table["n_orfs_annotated"] = table["n_orfs_annotated"].fillna(0).astype(int)

    in_neo = table["neordrp_best_orf"].notna()
    in_ips = table["ips_rdrp_best_orf"].notna()
    table.insert(1, "has_RdRp", (in_neo | in_ips).map({True: "yes", False: "no"}))
    evidence = pd.Series("-", index=table.index)
    evidence[in_neo & in_ips] = "both"
    evidence[in_neo & ~in_ips] = "NeoRdRp2"
    evidence[~in_neo & in_ips] = "InterProScan"
    table.insert(2, "rdrp_evidence", evidence)

    log(f"RdRp-positive contigs: {(table['has_RdRp'] == 'yes').sum()} "
        f"(both: {(evidence == 'both').sum()}, NeoRdRp2 only: {(evidence == 'NeoRdRp2').sum()}, "
        f"InterProScan only: {(evidence == 'InterProScan').sum()})")

    table.to_csv(args.out_table, sep="\t", index=False, na_rep="NA")

    # Best RdRp ORF per contig: NeoRdRp2 hit if present, otherwise the InterProScan one
    best_orf = table["neordrp_best_orf"].fillna(table["ips_rdrp_best_orf"]).dropna()
    with open(args.out_rdrp_faa, "w") as out:
        for orf_id in best_orf:
            out.write(f">{orf_id}\n{orfs[orf_id]}\n")
    log(f"Wrote {len(best_orf)} RdRp proteins to {args.out_rdrp_faa}")


if __name__ == "__main__":
    main()
