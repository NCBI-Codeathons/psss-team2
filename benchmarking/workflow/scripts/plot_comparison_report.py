import argparse
import json
import logging
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def check_ax(ax):
    if ax is None:
        return plt.gca()
    return ax


def plot_pr(data, ax=None, label=None, color=None):
    ax = check_ax(ax)
    ax.plot(data['qual_metrics']['recall'], data['qual_metrics']['precision'],
            label=label, color=color)
    ax.tick_params('both', labelsize='16')
    ax.set_xlabel('Recall', fontsize=18)
    ax.set_ylabel('Precision', fontsize=18)


def plot_recall(data, ax=None, label=None, color=None):
    ax = check_ax(ax)
    ax.plot(data['recall_pid']['pid'], data['recall_pid']['recall'],
            label=label, color=color)
    ax.tick_params('both', labelsize='16')
    ax.set_xlabel('Minmium average nucleotide identity', fontsize=18)
    ax.set_ylabel('Recall', fontsize=18)


def main(args):

    labels = args.labels
    if labels is None:
        labels = args.performance_reports
    else:
        labels = labels.split(',')

    if len(labels) != len(args.performance_reports):
        print((f'Please provide the same number of labels as reports. '
               f'Got {len(labels)} expected {len(args.performance_reports)}.'),
              file=sys.stderr)
        exit(1)

    if len(labels) > 20:
        print(f'You can plot at most 20 curves currently.')
        exit(1)

    pr_fig = plt.figure(figsize=(7, 7))
    rec_fig = plt.figure(figsize=(7, 7))

    cmap = mpl.cm.get_cmap('tab20')
    a = np.linspace(0, 1, 20)
    c = np.concatenate([a[0::2], a[1::2]])
    colors = [cmap(x) for x in c]

    for i, (report, label) in enumerate(zip(args.performance_reports, labels)):
        with open(report, 'r') as f:
            data = json.load(f)

        plot_pr(data, ax=pr_fig.gca(), label=label, color=colors[i])
        plot_recall(data, ax=rec_fig.gca(), label=label, color=colors[i])

    for fig, path in zip((pr_fig, rec_fig), ('pr_curve', 'recall_pid')):
        fig.tight_layout()
        fig.gca().legend()
        fig.savefig(f'{args.outdir}/{path}.png', dpi=200)


def parse_argparse_args():
    desc = """
    Generate plots for performance report
    """
    epi = """
    performance_report should be the output of performance_report.py
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("performance_reports", type=str, nargs='+',
                        help="performance reports in JSON format")
    parser.add_argument("-l", "--labels", type=str,
                        help="labels in comma-separated list", default=None)
    parser.add_argument("-o", "--outdir", type=str,
                        help="directory to save figures to", default='.')
    args = parser.parse_args()
    return args


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args.true_tsv = snakemake.input['performance_reports']
    args.outdir = snakemake.output.get('output_directory')
    args.labels = snakemake.arguments.get('labels')
    return args


if 'snakemake' in locals():
    args = parse_snakemake_args(snakemake)
    logging.basicConfig(
    filename=str(snakemake.log),
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
)
    logging.info(f"Starting script {__file__.split('/')[-1]}.")
    logging.debug(f"Full script path: {__file__}")
    main(args)
    logging.info(f"Done.")
elif __name__ == '__main__':
    args = parse_argparse_args()
    main(args)
