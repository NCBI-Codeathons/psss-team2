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
        get_novel_implementation_output(),
        get_performance_report_output(),
    )

    return outputs


def get_query_paths():
    return str(Path("output/").joinpath(config["query_paths"]))
    

def get_reference_paths():
    return str(Path("output/").joinpath(config["reference_paths"]))


def get_novel_implementation_output():
    return str(Path("output/").joinpath(config.get(
        "novel_implementation_output", "novel_implementation_output.tsv"
    )))


def get_performance_report_output():
    return str(Path("output/").joinpath(config.get(
        "performance_report_output", "performance_report_output.json"
    )))
