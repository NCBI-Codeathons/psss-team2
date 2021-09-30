from pathlib import Path


rule create_data_paths:
    input:
        data_table=config["data_table"],
    output:
        query_paths=get_query_paths(),
        reference_paths=get_reference_paths(),
    params:
        query=config["query"],
    log:
        "output/logs/create_data_paths.log",
    benchmark:
        "output/benchmarks/create_data_paths.txt"
    conda:
        "../envs/performance_report.yaml"
    script:
        "../scripts/create_data_paths.py"


rule download_s3_data:
    input:
        query_paths=get_query_paths(),
        reference_paths=get_reference_paths(),
    output:
        query_fna=get_query_fna(),
        reference_fna=get_reference_fna(),
    params:
        query_dir=lambda w, output: str(Path(output.query_fna).parent),
        reference_dir=lambda w, output: str(Path(output.reference_fna).parent),
    log:
        "output/logs/download_s3_data.log",
    benchmark:
        "output/benchmarks/download_s3_data.txt"
    conda:
        "../envs/aws.yaml"
    shell:
        """
        mkdir -p {params.query_dir} &>> {log}

        while read query_path; do
            aws s3 cp $query_path {params.query_dir} &>> {log}
        done < {input.query_paths}

        cat {params.query_dir}/*.fna > {output.query_fna} &>> {log}

        mkdir -p {params.reference_dir} &>> {log}

        while read reference_path; do
            aws s3 cp $reference_path {params.reference_dir} &>> {log}
        done < {input.reference_paths}

        cat {params.reference_dir}/*.fna > {output.reference_fna} &>> {log}
        """


rule bbmap_reformat:
    input:
        query_fna=get_query_fna(),
        reference_fna=get_reference_fna(),
    output:
        query_fna_filtered=get_query_fna_filtered(),
        reference_fna_filtered=get_reference_fna_filtered(),
    params:
        minlength=config["minlength"],
    log:
        "output/logs/bbmap_reformat.log",
    benchmark:
        "output/benchmarks/bbmap_reformat.txt"
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        reformat.sh in={input.query_fna} out={output.query_fna} ml={params.minlength} &>> {log}
        reformat.sh in={input.reference_fna} out={output.reference_fna} ml={params.minlength} &>> {log}
        """


rule samtools_faidx:
    input:
        query_fna_filtered=get_query_fna_filtered(),
        reference_fna_filtered=get_reference_fna_filtered(),
    output:
        query_fai=str(Path(get_query_fna_filtered()).with_suffix(".fai")),
        reference_fai=str(Path(get_reference_fna_filtered()).with_suffix(".fai")),
    threads: workflow.cores,
    log:
        "output/logs/samtools_faidx.log",
    benchmark:
        "output/benchmarks/samtools_faidx.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx -@ {threads} {input.query_fna_filtered} {output.query_fai} &>> {log}
        samtools faidx -@ {threads} {input.reference_fna_filtered} {output.reference_fai} &>> {log}
        """


rule mmseqs2:
    input:
        query_fna_filtered=get_query_fna_filtered(),
        reference_fna_filtered=get_reference_fna_filtered(),
    output:
        outfile_gz="output/mmseqs2_results.b6.gz",
    params:
        outfile=lambda w, output: output.outfile_gz[:-3],
        method="easy-search",
        search_type=3,
        tmp_dir="tmp"
    threads:
        workflow.cores,
    log:
        "output/logs/mmseqs2.log",
    benchmark:
        "output/benchmarks/mmseqs2.txt"
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs  {params.method}                     \
                --threads {threads}                 \
                --search-type {params.search_type}  \
                {input.query_fna_filtered}          \
                {input.reference_fna_filtered}      \
                {params.outfile}                    \
                {params.tmp_dir}

        gzip {params.outfile}
        """


rule novel_implementation:
    """
    Runs the novel implementation of the contig containment algorithm.

    Input:
        query -> contigs (in .FASTA format?),
        subject -> reference database 
    """
    input:
        query=get_query_fna_filtered(),
        subject=get_reference_fna_filtered(),
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


rule performance_report:
    input:
        predicted_containments="example/example_submission.tsv",
        ground_truth=config["ground_truth"],
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
