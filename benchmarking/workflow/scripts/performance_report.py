import argparse
import json
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sps


desc = """
Assess performance of a contig containment tool
"""

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("true_tsv", type=str,
                    help="the true containments in tabular BLAST output")
parser.add_argument("pred_tsv", type=str,
                    help="the predicted containments in tabular BLAST output")
parser.add_argument("-o", "--output", type=str,
                    help="the file to save results to", default=None)

args = parser.parse_args()

"""  Example TSV
qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore
nmdc:mga04781_15	nmdc:mga04781_2	97.6	8564	*	*	1	8565	5736	14300	*	*
nmdc:mga04781_3	nmdc:mga04781_15	95.8	6551	*	*	8865	15416	1	6552	*	*
"""

# read in files
true_df = pd.read_csv(args.true_tsv, sep='\t')
pred_df = pd.read_csv(args.pred_tsv, sep='\t')

# map files to indices
ctgs = set(true_df['qseqid'])
ctgs.update(true_df['sseqid'])

# find and filter extra contigs in predictions file
extras = set(pred_df['qseqid'])
extras.update(pred_df['sseqid'])
extras -= ctgs

if len(extras):
    print((f'Found {len(extras)} extra contigs in {args.pred_tsv}. '
            'Discarding before computing metrics.'), file=sys.stderr)

    mask = np.logical_or(pred_df['qseqid'].isin(ctgs),
                         pred_df['sseqid'].isin(ctgs))
    pred_df = pred_df[mask]

# build graph
## map sequence identifiers to indices
ctgs = dict(zip(ctgs, range(len(ctgs))))
true_df['qseqid_idx'] = [ctgs[s] for s in true_df['qseqid']]
true_df['sseqid_idx'] = [ctgs[s] for s in true_df['sseqid']]
pred_df['qseqid_idx'] = [ctgs[s] for s in pred_df['qseqid']]
pred_df['sseqid_idx'] = [ctgs[s] for s in pred_df['sseqid']]

## build adjacency matrix
n_ctgs = len(ctgs)
true_g = sps.dok_matrix((n_ctgs, n_ctgs), dtype=np.int8)
pred_g = sps.dok_matrix((n_ctgs, n_ctgs), dtype=np.int8)
true_g[true_df['qseqid_idx'], true_df['sseqid_idx']] = 1
pred_g[pred_df['qseqid_idx'], pred_df['sseqid_idx']] = 1
true_g = true_g.tocsr()
pred_g = pred_g.tocsr()

# calculate metrics:
diff = true_g - pred_g
fp = (diff < 0).sum()
fn = (diff > 0).sum()
tp = true_g.sum() - fn
tn = n_ctgs**2 - tp - fn - fp

metrics = dict()
metrics['precision'] = tp / (tp + fp)
metrics['sensitivity'] = tp / (tp + fn)
metrics['specificity'] = tn / (tn + fp)
metrics['accuracy'] = (tp + tn) / (tn + tp + fn + fp)
metrics['f1'] = tp / (tp + (fp + fn)/2)

# output metrics
if args.output is not None:
    out = open(args.output, 'w')
else:
    out = sys.stdout

json.dump(metrics, fp=out)