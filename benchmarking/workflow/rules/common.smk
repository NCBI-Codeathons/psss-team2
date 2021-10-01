"""
common.smk - helper functions for other rules.

Contrary to the other .smk files, this does not use the Snakemake DSL and rather it's just usual
Python syntax.
"""

from pathlib import Path


# Validators
def validate_inputs():

    # These must exist in the config.yaml file
    inputs = ("subject", "ground_truth")

    for k in inputs:
        file = config.get(k, "")
        assert Path(file).exists(), f"Input file '{k}' not found: {file}."


# Validator calls
validate_inputs()


# Outputs
def get_final_output():
    outputs = (
        "output/mmseqs2_results.b6.gz",
        "output/blastn_results.b6.gz",
        "output/synteny_results.b6.gz",
        get_novel_implementation_output(),
        get_performance_report_output(),
    )

    return outputs


def assert_valid_data_kind(kind):
    valid_kinds = ("query", "reference")
    assert kind in valid_kinds, f"'kind' must be one of {valid_kinds}, got '{kind}'."


def get_paths(kind="query"):
    assert_valid_data_kind(kind)
    return str(Path("output/").joinpath(config[f"{kind}_paths"]))


def get_fna(kind="query"):
    assert_valid_data_kind(kind)
    return str(Path(config["data_dir"]).joinpath(f"{kind}/{kind}.fna"))


def get_fna_filtered(kind="query"):
    assert_valid_data_kind(kind)
    minlength = config.get("minlength", 500)
    return str(Path(get_fna(kind)).with_suffix(f".ml{minlength}.fna"))


def get_fai(kind="query"):
    assert_valid_data_kind(kind)
    return str(Path(get_fna_filtered(kind)).with_suffix(".fai"))

def get_makeblastdb_out():
    return [str(Path(get_fna_filtered("reference")).with_suffix(".fna.nhr")), ]


def get_novel_implementation_output():
    return str(
        Path("output/").joinpath(
            config.get("novel_implementation_output", "novel_implementation_output.tsv")
        )
    )


def get_performance_report_output():
    return str(
        Path("output/").joinpath(
            config.get("performance_report_output", "performance_report_output.json")
        )
    )
