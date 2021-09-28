
rule novel_implementation:
    """
    Runs the novel implementation of the contig containment algorithm.

    Input:
        query -> contigs (in .FASTA format?),
        subject -> reference database 
    """
    input:
        query=config["query"],
        subject=config["subject"],
    output:
        outfile=get_novel_implementation_output(),
    log:
        "output/logs/novel_implementation.log",
    benchmark:
        "output/benchmarks/novel_implementation.txt"
    conda:
        "../envs/novel_implementation.yaml"
    shell:
        # Could look like this: {{ time -v novel_implementation {input.query} {input.subject} ; }} &> {log}
        # Will only be mock for now
        """
        touch {output}
        """


rule post_processing:
    input:
        novel_implementation_outfile=get_novel_implementation_output(),
        ground_truth=config["ground_truth"],
    output:
        outfile=get_post_processing_output(),
    log:
        "output/logs/post_processing.log",
    benchmark:
        "output/benchmarks/post_processing.txt"
    conda:
        "../envs/post_processing.yaml"
    shell:
        """
        touch {output}
        """


rule performance_report:
    input:
        post_processing_outfile=get_post_processing_output(),
    output:
        performance_report=get_performance_report_output(),
    log:
        "output/logs/performance_report.log",
    benchmark:
        "output/benchmarks/performance_report.txt"
    conda:
        "../envs/performance_report.yaml"
    script:
        "../scripts/performance_report.py"
