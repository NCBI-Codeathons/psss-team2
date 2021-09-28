# Benchmarking code for identifying contig similarity across metagenomic datasets
Last modified on: September 28, 2021

## Introduction
This project contains code for evaluating the performance and accuracy 
of tools for identifying the similarity between one or more query contigs and 
a "reference" collection of metagenomic samples.

In the initial implementation, the focus is solely on containment relationships.

## Limitations
The task being benchmarked involves comparisons between **contigs** rather than
individual reads.  Thus, both the query and the reference are assumed to be 
derived from assembled samples and represented as DNA strings that are relatively
long (e.g., longer than 500 basepairs)

## Inputs
The benchmarking tool takes as input:
* A "gold standard" output in TAB delimited format as described in more detail below
* The file paths to the query and reference contigs sets. The reference may comprise one
  or more paths to metagenomic samples.
* The path to an executable that will run the code being benchmarked

These inputs are provided in [??Yaml file??]

## Outputs
The benchmarking tool runs the code to be benchmarked then produces a report that includes
several metrics describing the "match" between the gold standard and the output of the tool.

[If time permits] Information about each run is stored in a database in order to allow
comparisons across tools and/or data sets.

## Problem definition

 
## Expectations from benchmarked tool
The benchmarked tool must accept two parameters

## Metrics measured
