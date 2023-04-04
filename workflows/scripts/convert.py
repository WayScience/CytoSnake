"""
script: covert.py

converts sqlite (or other formats) into parquet files
"""
import pathlib
from typing import Union, List

import yaml
import cytotable


def sqlite_to_parquet(
    input_file: Union[str, List[str]],
    out_path: str,
    target_ext: str,
    convert_config_path: str,
):
    # loading config files and only selecting the parameters
    convert_config_path = pathlib.Path(convert_config_path).resolve(strict=True)
    with open(convert_config_path, "r") as yaml_contents:
        cytotable_config = yaml.safe_load(yaml_contents)["cytotable_convert"]["params"]

    # type checking, making sure that a sqlite file is passed.
    if target_ext == ".parquet" and pathlib.Path(input_file).suffix != ".sqlite":
        raise ValueError(
            "Converting to parquet files requires sqlite file."
            f"File provided: {target_ext}"
        )

    # convert sqlite file into parquet
    cytotable.convert(
        source_path=input_file,
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
    plate_data = str(snakemake.input)
    output_path = str(snakemake.output)
    general_configs = snakemake.params["data_configs"]
    config_path = snakemake.params["cytotable_config"]

    # executing conversion input file to parquet
    sqlite_to_parquet(
        input_file=plate_data,
        out_path=output_path,
        convert_config_path=config_path,
        target_ext=general_configs,
    )


if __name__ == "__main__":
    main()
