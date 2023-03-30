"""
script: covert.py

converts sqlite (or other formats) into parquet files
"""
import pathlib
from typing import Union

import cytotable

from cytosnake.guards.ext_guards import has_sqlite_ext
from cytosnake.utils.config_utils import load_configs, load_general_configs
from cytosnake.common.errors import ExtensionError


def sqlite_to_parquet(
    sqlite_file: Union[str, list[str]],
    out_path: str,
    config_file: Union[str, pathlib.Path],
):
    # loading config files and only selecting the parameters
    general_configs = load_general_configs()
    cytotable_config = load_configs(config_file)["params"]

    # type checking, making sure that a sqlite file is passed.
    target_ext = general_configs["data_configs"]["plate_data_format"]
    if target_ext == ".parquet" and not has_sqlite_ext(sqlite_file):
        raise ExtensionError(
            "Converting to parquet files requires sqlite file."
            f"File provided: {target_ext}"
        )

    # convert sqlite file into parquet
    cytotable.convert(
        source_path=sqlite_file,
        dest_path=out_path,
        dest_datatype=cytotable_config["dest_datatype"],
        source_datatype=cytotable_config["source_datatype"],
        concat=cytotable_config["concat"],
        join=cytotable_config["join"],
        infer_common_schema=cytotable_config["infer_common_schema"],
        drop_null=cytotable_config["drop_null"],
        preset=cytotable_config["preset"],
        log_level=cytotable_config["log_level"],
    )


def main():
    """Execution of the main script"""

    # grabbing snakemake inputs
    plate_data = snakemake.input
    output_path = snakemake.output
    config_path = snakemake.params["cytotable_config"]

    # executing conversion input file to parquet
    sqlite_to_parquet(
        source_path=plate_data, dest_path=output_path, config_path=config_path
    )


if __name__ == "__main__":
    main()
