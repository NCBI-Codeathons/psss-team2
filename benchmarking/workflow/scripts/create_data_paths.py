import pandas as pd


def make_S3_path(sample):
    return "s3://psss-metagenomics-codeathon-data/marine/%s/assembly/%s_contigs.fna" % (
        sample,
        sample.replace(":", "_"),
    )


def main(snakemake):
    QUERY_SAMPLE = snakemake.params["query"]
    
    metadata = pd.read_csv(snakemake.input["data_table"], sep="\t")
    metadata = metadata.query("TAXA == 'marine metagenome'")
    metadata = metadata.query("MGA_ID != @QUERY_SAMPLE")
    full_paths = [make_S3_path(i) for i in metadata["MGA_ID"].values]


    with open(snakemake.output["query_paths"], "w") as f:
        f.write(make_S3_path(QUERY_SAMPLE) + "\n")

    with open(snakemake.output["reference_paths"], "w") as f:
        f.write("%s\n" % "\n".join(full_paths))


if "snakemake" in locals() and __name__ == "__main__":
    main(snakemake)
