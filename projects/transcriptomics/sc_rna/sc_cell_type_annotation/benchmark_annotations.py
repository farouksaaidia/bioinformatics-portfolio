#!/usr/bin/env python3
"""
Benchmarking annotations:
- Inputs: gold CSV (cell_id,label_gold), pred CSV (cell_id,label_pred)
- Outputs: classification report CSV (precision/recall/f1 per label), normalized confusion matrix CSV, and a small PNG heatmap for misclassifications.
"""
import argparse
import os
import sys
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--gold", required=True, help="Gold labels CSV (cell_id,label_gold)")
parser.add_argument("--pred", required=True, help="Pred labels CSV (cell_id,label_pred)")
parser.add_argument("--out_prefix", required=True, help="Output file prefix")
args = parser.parse_args()

for f in [args.gold, args.pred]:
    if not os.path.exists(f):
        sys.exit(f"❌ {f} not found")

gold = pd.read_csv(args.gold)
pred = pd.read_csv(args.pred)
if not {'cell_id','label_gold'}.issubset(gold.columns) or not {'cell_id','label_pred'}.issubset(pred.columns):
    sys.exit("❌ gold must have (cell_id,label_gold), pred must have (cell_id,label_pred)")

merged = pd.merge(gold, pred, on='cell_id')
labels = sorted(list(set(merged.label_gold.unique()) | set(merged.label_pred.unique())))

prfs = precision_recall_fscore_support(merged.label_gold, merged.label_pred, labels=labels, zero_division=0)
report = pd.DataFrame({'label': labels, 'precision': prfs[0], 'recall': prfs[1], 'f1': prfs[2], 'support': prfs[3]})
report.to_csv(f"{args.out_prefix}_classification_report.csv", index=False)

cm = confusion_matrix(merged.label_gold, merged.label_pred, labels=labels)
cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
pd.DataFrame(cm, index=labels, columns=labels).to_csv(f"{args.out_prefix}_confusion.csv")
pd.DataFrame(cm_norm, index=labels, columns=labels).to_csv(f"{args.out_prefix}_confusion_normalized.csv")

# heatmap
plt.figure(figsize=(8,6))
plt.imshow(cm_norm, interpolation='nearest')
plt.title('Normalized Confusion Matrix')
plt.colorbar()
plt.xticks(range(len(labels)), labels, rotation=90)
plt.yticks(range(len(labels)), labels)
plt.tight_layout()
plt.savefig(f"{args.out_prefix}_confusion_normalized.png", dpi=150)
print(f"✅ Benchmark files written with prefix {args.out_prefix}")
