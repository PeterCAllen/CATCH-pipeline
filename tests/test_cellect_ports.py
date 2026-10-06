#!/usr/bin/env python3
"""Synthetic-data tests for the CELLECT ports and the specificity sanitizer."""
import gzip
import os
import subprocess
import sys
import tempfile

import pandas as pd

BIN = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin")
FAILED = []


def run(script, *args):
    cmd = [sys.executable, os.path.join(BIN, script)] + list(args)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(f"{script} failed:\nSTDOUT:{r.stdout}\nSTDERR:{r.stderr}")
    return r.stdout


def check(name, cond, detail=""):
    if cond:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILED.append(name)


def test_format_genes(d):
    print("\n[cellect_format_genes]")
    coords = os.path.join(d, "coords.txt")
    with open(coords, "w") as fh:
        # GENE CHR START END STRAND GENE_NAME
        fh.write("ENSG1\t1\t500\t1500\t+\tA\n")        # start clipped to 0
        fh.write("ENSG2\t1\t400000\t401000\t-\tB\n")
        fh.write("ENSG3\t2\t250\t900\t+\tC\n")         # END clipped to chr size
        fh.write("ENSG4\t23\t100\t200\t+\tX1\n")       # chr not in sizes -> dropped

    sizes = os.path.join(d, "sizes.txt")
    with open(sizes, "w") as fh:
        fh.write("1\t1000000\n2\t100500\n")

    out = os.path.join(d, "genes")
    os.makedirs(out, exist_ok=True)
    run("cellect_format_genes.py", "--gene_coords", coords, "--chr_sizes", sizes,
        "--out_dir", out, "--windowsize_kb", "100")

    b1 = pd.read_csv(os.path.join(out, "genes_plus_100kb.1.bed"), sep="\t", header=None,
                     names=["CHR", "START", "END", "GENE"])
    b2 = pd.read_csv(os.path.join(out, "genes_plus_100kb.2.bed"), sep="\t", header=None,
                     names=["CHR", "START", "END", "GENE"])

    check("chr1 has 2 genes", len(b1) == 2, f"got {len(b1)}")
    check("START floors at 0", int(b1.loc[b1.GENE == "ENSG1", "START"].iloc[0]) == 0)
    check("END = end + 100kb", int(b1.loc[b1.GENE == "ENSG1", "END"].iloc[0]) == 101500)
    check("ENSG2 START = 300000", int(b1.loc[b1.GENE == "ENSG2", "START"].iloc[0]) == 300000)
    check("END clips to chr size", int(b2.loc[b2.GENE == "ENSG3", "END"].iloc[0]) == 100500)
    check("unknown chr dropped", not os.path.exists(os.path.join(out, "genes_plus_100kb.23.bed")))
    check("sorted by START", list(b1.START) == sorted(b1.START))
    check("headerless BED4", b1.shape[1] == 4)


def test_make_all_genes(d):
    print("\n[cellect_make_all_genes]")
    spec = os.path.join(d, "spec.csv")
    pd.DataFrame({"gene": ["G1", "G2", "G3"], "ct_a": [0.1, 0.0, 0.9],
                  "ct_b": [0.5, 0.2, 0.0]}).to_csv(spec, index=False)
    out = os.path.join(d, "all_genes.csv")
    run("cellect_make_all_genes.py", "--spec_matrix", spec, "--out", out)

    got = pd.read_csv(out)
    check("has gene + control column", list(got.columns) == ["gene", "all_genes_in_dataset"],
          str(list(got.columns)))
    check("all values are 1", set(got.all_genes_in_dataset) == {1})
    check("one row per gene", len(got) == 3)


def test_make_cts_file(d):
    print("\n[cellect_make_cts_file]")
    ann = os.path.join(d, "annots.txt")
    with open(ann, "w") as fh:
        fh.write("ct_a\nct_b\n\n")
    out = os.path.join(d, "x.ldcts")
    run("cellect_make_cts_file.py", "--annotations", ann, "--run_prefix", "cepo_norm",
        "--ldscore_dir", ".", "--out", out)

    lines = [l.rstrip("\n") for l in open(out) if l.strip()]
    check("one line per annotation", len(lines) == 2, str(lines))
    name, prefix = lines[0].split("\t")
    check("name is <prefix>__<annot>", name == "cepo_norm__ct_a", name)
    check("path ends with a dot", prefix.endswith("."), prefix)
    check("blank annotation lines ignored", all(l.count("\t") == 1 for l in lines))


def test_parse_results(d):
    print("\n[cellect_parse_results]")
    ctr = os.path.join(d, "res.cell_type_results.txt")
    pd.DataFrame({
        "Name": ["cepo_norm__ct_a", "cepo_norm__ct_b"],
        "Coefficient": [1.2e-8, -3.0e-9],
        "Coefficient_std_error": [4.0e-9, 5.0e-9],
        "Coefficient_P_value": [0.0013, 0.72],
    }).to_csv(ctr, sep="\t", index=False)

    out = os.path.join(d, "prioritization.csv")
    run("cellect_parse_results.py", "--cell_type_results", ctr, "--gwas", "TEST", "--out", out)

    got = pd.read_csv(out)
    check("CELLECT schema",
          list(got.columns) == ["gwas", "specificity_id", "annotation", "beta", "beta_se", "pvalue"],
          str(list(got.columns)))
    check("specificity_id split", set(got.specificity_id) == {"cepo_norm"})
    check("annotation split", set(got.annotation) == {"ct_a", "ct_b"})
    check("pvalue passthrough", abs(got.loc[got.annotation == "ct_a", "pvalue"].iloc[0] - 0.0013) < 1e-12)
    check("sorted by pvalue", list(got.pvalue) == sorted(got.pvalue))


def test_split_ldscores(d):
    print("\n[cellect_split_ldscores]")
    base = os.path.join(d, "cepo_norm.COMBINED_ANNOT.7")
    ld = pd.DataFrame({
        "CHR": [7, 7, 7], "SNP": ["rs1", "rs2", "rs3"], "BP": [10, 20, 30],
        "ct_aL2": [0.5, 0.6, 0.7], "ct_bL2": [1.5, 1.6, 1.7],
    })
    ld.to_csv(base + ".l2.ldscore.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame({"ct_a": [11.0], "ct_b": [22.0]}).to_csv(
        base + ".annot.gz", sep="\t", index=False, compression="gzip")
    with open(base + ".l2.M", "w") as fh:
        fh.write("11.0\t22.0\n")
    with open(base + ".l2.M_5_50", "w") as fh:
        fh.write("9.0\t18.0\n")

    outdir = os.path.join(d, "per_annotation")
    run("cellect_split_ldscores.py", "--ldscore", base + ".l2.ldscore.gz", "--out_dir", outdir)

    fa = os.path.join(outdir, "cepo_norm__ct_a.7.l2.ldscore.gz")
    check("per-annotation ldscore written", os.path.exists(fa))
    sub = pd.read_csv(fa, sep="\t")
    check("ldscore keeps CHR/SNP/BP + one annot",
          list(sub.columns) == ["CHR", "SNP", "BP", "ct_aL2"], str(list(sub.columns)))
    check("M split to correct field",
          open(os.path.join(outdir, "cepo_norm__ct_b.7.l2.M")).read().strip() == "22.0")
    check("M_5_50 split to correct field",
          open(os.path.join(outdir, "cepo_norm__ct_a.7.l2.M_5_50")).read().strip() == "9.0")
    with gzip.open(os.path.join(outdir, "cepo_norm__ct_a.7.annot.gz"), "rt") as fh:
        hdr = fh.readline().strip()
    check("annot split keeps header", hdr == "ct_a", hdr)


def test_sanitize(d):
    print("\n[sanitize_specificity_matrix]")
    src = os.path.join(d, "raw.csv")
    pd.DataFrame({
        "Unnamed: 0": ["G1", "G2", "G3", "G4"],
        "ct_a": [0.9, None, 0.1, 0.4],
        "ct_zero": [0.0, 0.0, 0.0, 0.0],
        "ct_b": [0.2, 0.8, 0.3, 0.5],
    }).to_csv(src, index=False)

    out = os.path.join(d, "clean.csv")
    ann = os.path.join(d, "annots.txt")
    run("sanitize_specificity_matrix.py", "--input", src, "--output", out,
        "--output_annotations", ann)

    got = pd.read_csv(out)
    check("first column renamed to gene", got.columns[0] == "gene", got.columns[0])
    check("all-zero annotation dropped", "ct_zero" not in got.columns)
    check("NA filled with 0", float(got.loc[got.gene == "G2", "ct_a"].iloc[0]) == 0.0)
    check("annotation list matches", [l.strip() for l in open(ann) if l.strip()] == ["ct_a", "ct_b"])

    out2 = os.path.join(d, "clean_top.csv")
    run("sanitize_specificity_matrix.py", "--input", src, "--output", out2,
        "--output_annotations", os.path.join(d, "a2.txt"), "--top_pct", "0.5")
    got2 = pd.read_csv(out2)
    check("top_pct zeroes the bottom half",
          int((got2.ct_b > 0).sum()) == 2, f"{list(got2.ct_b)}")
    check("top_pct keeps the largest values",
          set(got2.loc[got2.ct_b > 0, "gene"]) == {"G2", "G4"}, str(list(got2.gene[got2.ct_b > 0])))

    bad = os.path.join(d, "bad.csv")
    pd.DataFrame({"gene": ["G1"], "bad__name": [0.5]}).to_csv(bad, index=False)
    r = subprocess.run([sys.executable, os.path.join(BIN, "sanitize_specificity_matrix.py"),
                        "--input", bad, "--output", os.path.join(d, "o.csv"),
                        "--output_annotations", os.path.join(d, "a.txt")],
                       capture_output=True, text=True)
    check("rejects '__' in annotation names", r.returncode != 0 and "CELLECT-safe" in r.stderr)


with tempfile.TemporaryDirectory() as d:
    test_format_genes(d)
    test_make_all_genes(d)
    test_make_cts_file(d)
    test_parse_results(d)
    test_split_ldscores(d)
    test_sanitize(d)

print()
if FAILED:
    print(f"{len(FAILED)} FAILED: {FAILED}")
    sys.exit(1)
print("All tests passed.")
