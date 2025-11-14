#!/usr/bin/env python3
import argparse, os, sys
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np

parser = argparse.ArgumentParser(description="Identify short-read RNA-seq outlier samples via PCA + MAD")
parser.add_argument("-c","--counts", required=True, help="Normalized counts file")
parser.add_argument("-o","--output", required=True, help="Output TSV with outlier flags")
args = parser.parse_args()

if not os.path.exists(args.counts):
    sys.exit("Counts file not found.")

df = pd.read_csv(args.counts, sep=",", index_col=0)
mat = df.values.T

scaled = StandardScaler().fit_transform(mat)
pca = PCA(n_components=2).fit_transform(scaled)

dists = np.sqrt((pca[:,0]**2) + (pca[:,1]**2))
mad = np.median(np.abs(dists - np.median(dists)))
threshold = np.median(dists) + 4 * mad

outliers = dists > threshold

result = pd.DataFrame({
    "sample": df.columns,
    "pca_distance": dists,
    "outlier": outliers
})

result.to_csv(args.output, sep="\t", index=False)
print("Outlier detection complete.")
