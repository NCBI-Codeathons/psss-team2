from pathlib import Path


rule create_data_paths:
    input:
        data_table=config["data_table"],
    output:
        query_paths=get_paths("query"),
        reference_paths=get_paths("reference"),
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
        query_paths=get_paths("query"),
        reference_paths=get_paths("reference"),
    output:
        query_fna=get_fna("query"),
        reference_fna=get_fna("reference"),
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

        {{ cat {params.query_dir}/*.fna > {output.query_fna} ; }} &>> {log}

        mkdir -p {params.reference_dir} &>> {log}

        while read reference_path; do
            aws s3 cp $reference_path {params.reference_dir} &>> {log}
        done < {input.reference_paths}

        {{ cat {params.reference_dir}/*.fna > {output.reference_fna} ; }} &>> {log}
        """


rule bbmap_reformat:
    input:
        query_fna=get_fna("query"),
        reference_fna=get_fna("reference"),
    output:
        query_fna_filtered=get_fna_filtered("query"),
        reference_fna_filtered=get_fna_filtered("reference"),
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
        reformat.sh in={input.query_fna} out={output.query_fna_filtered} ml={params.minlength} &>> {log}
        reformat.sh in={input.reference_fna} out={output.reference_fna_filtered} ml={params.minlength} &>> {log}
        """


rule samtools_faidx:
    input:
        query_fna_filtered=get_fna_filtered("query"),
        reference_fna_filtered=get_fna_filtered("reference"),
    output:
        query_fai=get_fai("query"),
        reference_fai=get_fai("reference"),
    log:
        "output/logs/samtools_faidx.log",
    benchmark:
        "output/benchmarks/samtools_faidx.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx {input.query_fna_filtered} &>> {log}
        samtools faidx {input.reference_fna_filtered} &>> {log}
        """


rule mmseqs2:
    input:
        query_fna_filtered=get_fna_filtered("query"),
        query_fai=get_fai("query"),
        reference_fna_filtered=get_fna_filtered("reference"),
        reference_fai=get_fai("reference"),
    output:
        outfile_gz="output/mmseqs2_results.b6.gz",
        outfile_filtered_gz="output/mmseqs2_results_filtered.b6.gz",
    params:
        outfile=lambda w, output: output.outfile_gz[:-3],
        outfile_filtered=lambda w, output: output.outfile_filtered_gz[:-3],
        method="easy-search",
        search_type=3,
    threads: workflow.cores
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
                {params.outfile} /tmp &>> {log}

        python workflow/scripts/filter_blast_6_to_containments.py   \
               -i {params.outfile}                                  \
               -q {input.query_fai}                                 \
               -r {input.reference_fai}                             \
               -o {params.outfile_filtered}                         \
               --as_decimal &>> {log}

        gzip {params.outfile} &>> {log}
        gzip {params.outfile_filtered} &>> {log}
        """


rule synteny:
    input:
        query_fna_filtered=get_fna_filtered("query"),
        reference_fna_filtered=get_fna_filtered("reference"),
    output:
        synteny_outfile_gz="output/synteny_results.b6.gz",
    params:
        synteny_outfile=lambda w, output: output.synteny_outfile_gz[:-3],
        reference_folder=lambda w, input: str(
            Path(input.reference_fna_filtered).parent
        ),
    log:
        "output/logs/synteny.log",
    benchmark:
        "output/benchmarks/synteny.txt"
    conda:
        "../envs/synteny.yaml"
    shell:
        """
        Rscript workflow/scripts/run_synteny.R          \
                {params.reference_folder}               \
                {input.query_fna_filtered}              \
                {params.synteny_outfile} &> {log}

        gzip {params.synteny_outfile} &>> {log}
        """


rule makeblastdb:
    input:
        reference_fna_filtered=get_fna_filtered("reference"),
    output:
        outfiles=get_makeblastdb_out(),
    params:
        input_type="fasta",
        dbtype="nucl",
    log:
        "output/logs/makeblastdb.log",
    benchmark:
        "output/benchmarks/makeblastdb.txt"
    conda:
        "../envs/blast.yaml"
    shell:
        """
        makeblastdb -in {input}                     \
                    -input_type {params.input_type} \
                    -dbtype {params.dbtype}         \
                    -parse_seqids &> {log}
        """


rule novel_implementation:
    """
    Runs the novel implementation of the contig containment algorithm.

    Input:
        query -> contigs (in .FASTA format?),
        subject -> reference database 
    """
    input:
        query=get_fna_filtered("query"),
        subject=get_fna_filtered("reference"),
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


rule blastn:
    input:
        query=get_fna_filtered("query"),
        query_fai=get_fai("query"),
        reference_fai=get_fai("reference"),
        makeblastdb_out=get_makeblastdb_out(),
    output:
        outfile_gz="output/blastn_results.b6.gz",
        outfile_filtered_gz="output/blastn_results_filtered.b6.gz",
    params:
        outfmt=6,
        outfile=lambda w, output: output.outfile_gz[:-3],
        outfile_filtered=lambda w, output: output.outfile_filtered_gz[:-3],
        db=lambda w, input: str(Path(input.makeblastdb_out[0]).with_suffix("")),
    threads: workflow.cores
    log:
        "output/logs/blastn.log",
    benchmark:
        "output/benchmarks/blastn.txt"
    conda:
        "../envs/blast.yaml"
    shell:
        """
        blastn  -db {params.db}                 \
                -query {input.query}            \
                -out {params.outfile}           \
                -outfmt {params.outfmt}         \
                -num_threads {threads}          \
                &>> {log}

        python workflow/scripts/filter_blast_6_to_containments.py   \
               -i {params.outfile}                                  \
               -q {input.query_fai}                                 \
               -r {input.reference_fai}                             \
               -o {params.outfile_filtered} &>> {log}

        gzip {params.outfile} &>> {log}
        gzip {params.outfile_filtered} &>> {log}
        """


rule performance_report:
    input:
        predicted_containments="output/mmseqs2_results_filtered.b6.gz",
        ground_truth="output/blastn_results_filtered.b6.gz",
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
