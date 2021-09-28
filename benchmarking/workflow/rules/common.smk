"""
common.smk - helper functions for other rules.

Contrary to the other .smk files, this does not use the Snakemake DSL and rather it's just usual
Python syntax.
"""

from pathlib import Path


# Validators
def validate_inputs():

    # These must exist in the config.yaml file
    inputs = (
        "query",
        "subject",
        "ground_truth"
    )

    for k in inputs:
        file = config.get(k)
        assert Path(file).exists(), f"Input file '{k}' not found: {file}."


# Validator calls
validate_inputs()


# Outputs
def get_final_output():
    outputs = (
        get_novel_implementation_output(),
        get_post_processing_output(),
        get_performance_report_output()
    )

    return output


def get_novel_implementation_output():
    return config.get("novel_implementation")

def get_post_processing_output():
    return config.get("post_processing")

def get_performance_report_output():
    return config.get("performance_report")

