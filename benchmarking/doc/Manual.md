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
For the purposes of this application, we define a "match" between two contigs as one of the
following cases:

### Containment
Given contigs A of length lenA and B of length lenB, we say that A and B have a containment
relationship if an alignment exists between A and B that has higher than 95% identity within the aligned
region, and that covers more than 95% of the length of the shortest contig. 

A containment should be reported irrespective of whether the query is the longer or the shorter
among the contigs in the containment.  

Assume the length of contig A is 1000 bp and the length of contig B is 10000, the following
are examples of valid containment relationships.

[TODO - include several examples of tab-delimited output for correct and incorrect matches]
[TODO - all of A included]
[TODO - < 5% of A mismatches at one or both ends]

Here are several examples of relationships that don't represent valid containments
[TODO - alignment < 95% ]
[TODO - < 95% of A matches with mismatches at one or both ends]

### Overlap
[TODO] to be defined later

### Inner match [TODO - get better name]
[TODO] to be defined later
 
## Expectations from benchmarked tool
The benchmarked tool must accept two parameters: a path to the query metagenome and
a path to one or more reference samples. 

## Tab-delimited representation of matches
We are using a modified version of the BLAST tabular output (`--outfmt 6`). BLAST 
outputs the following fields, separated by `TAB` characters:

   1. **qseqid**:
      query sequence id

   2. **sseqid**:      
      reference sequence id

   3. **pident**:
      percentage of identical matches within the aligned region

   4. **length**:      
      alignment length

   5. **mismatch**:    
      number of mismatches

   6. **gapopen**:     
      number of gap openings

   7. **qstart**:
      start of alignment in query (1-based)

   8. **qend**:
       end of alignment in query (1-based/inclusive)

   9. **sstart**:      
       start of alignment in reference sequence (1-based)

   10. **send**:        
       end of alignment in subject (1-based/inclusive)

   11. **evalue**:      
       expect value

   12. **bitscore**:    
       bit score

The gold standard output must include field 1-4,7-10, with the other fields optional.
The benchmarked tool must include fields 1 and 2 with the other fields optional (i.e.,
it suffices to report a relationship exists between the contigs).

The files must only include the rows that correspond to valid matches between the 
contigs.

## Metrics measured
The following metrics are computed:

**true positive (TP)**:
  row in the output of the benchmarked tool that matches a roww in the gold
standard file

**true negative (TN)**:
  pair of contigs that is missing from the output of both the gold standard
and benchmarked tool

**false positive (FP)**:
  row in the output of the benchmarked tool that does not have a corresponding
row in the gold standard

**false negative (FN)**:
  row in the gold standard that does not have a corresponding row in the output
of the benchmarked tool

**precision**:
  $$precision = \frac{TP}{TP+FP}$$

**recall**:
  $$recall = \frac{TP}{TP+FN}$$

**maximum memory usage**:
  [TODO] definitions for single thread, multi-thread, and distributed computation versions

**runtime**:
  [TODO] definitions for single thread, multi-thread, and distributed computation versions

Additional reports stratify the accuracy measures described above by the percent identity
in the gold standard and by the length of the smaller contig in each pair within the
gold standard.
