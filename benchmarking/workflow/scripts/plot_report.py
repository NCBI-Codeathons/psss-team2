import argparse
import json
import logging
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def main(args):

    with open(args.performance_report, 'r') as f:
        data = json.load(f)

    p_recall = data['recall_pid']['recall']
    pid = data['recall_pid']['pid']

    plt.figure(figsize=(7, 7))
    ax = plt.gca()
    ax.plot(pid, p_recall, c='k')
    ax.tick_params('both', labelsize='16')
    ax.set_xlabel('Average Nucleotide Identity', fontsize=18)
    ax.set_ylabel('Recall', fontsize=18)
    plt.savefig(f'{args.outdir}/recall_by_pid.png')

    plt.figure(figsize=(7, 7))
    ax = plt.gca()
    q_recall = data['qual_metrics']['recall']
    q_precision = data['qual_metrics']['precision']
    min_qual = data['qual_metrics']['min_qual']

    ax.plot(pid, q_recall, c='k', label='Recall')
    ax.plot(pid, q_precision, c='r', label='Precision')
    ax.tick_params('both', labelsize='16')
    ax.set_xlabel('Prediction Quality', fontsize=18)
    ax.legend(loc='lower left', fontsize=16)
    plt.savefig(f'{args.outdir}/metrics_by_quality.png')


def parse_argparse_args():
    desc = """
    Generate plots for performance report
    """
    epi = """
    performance_report should be the output of performance_report.py
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("performance_report", type=str,
                        help="performance report in JSON format")
    parser.add_argument("-o", "--outdir", type=str,
                        help="directory to save figures to", default='.')
    args = parser.parse_args()
    return args


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args.true_tsv = snakemake.input['performance_report']
    args.outdir = snakemake.output['output_directory']
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
