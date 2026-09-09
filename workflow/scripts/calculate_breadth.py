#!/usr/bin/env python3
import pandas as pd
import sys

cov_file = sys.argv[1]
out_file = sys.argv[2]

df = pd.read_csv(cov_file, sep="\t", header=None,
                 names=["contig", "pos", "depth"])

covered = df[df["depth"] > 0].groupby("contig").size()
lengths = df.groupby("contig")["pos"].max()

breadth = (covered / lengths).fillna(0)

breadth.to_csv(out_file, sep="\t", header=True)
