"""
script: convert.py

converts sqlite (or other formats) into parquet files
"""
import pathlib
from typing import List, Union

import cytotable


def convert_to_parquet(
    input_file: Union[str, List[str]],
    out_path: str,
    target_ext: str,
    convert_configs: str,
):
    """Takes in a file or a list of files that will be converted into parquet.

    Parameters
    ----------
    input_file : Union[str, List[str]]
        files or list of files converted into `.parquet` files
    out_path : str
        path where generated parquet files will be saved
    target_ext : str
        dictates which file format
    convert_configs : dict
        dictionary containing cytotable.convert() configs

    Raises
    ------
    ValueError
        raised if an unsupported file extension is provided
    """
    # checking if user has parquet file
    if target_ext == ".parquet" and pathlib.Path(input_file).suffix not in [
        ".sqlite",
        ".csv",
    ]:
        raise ValueError(
            "Converting to parquet files requires sqlite file."
            f"File provided: {target_ext}"
        )

    # convert sqlite file into parquet
    cytotable_config = convert_configs["params"]
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

    # grabbing snakemake inputs from workflow
    # more information how snakemake transfers workflow variables to scripts:
    # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#python
    plate_data = str(snakemake.input)
    output_path = str(snakemake.output)
    general_configs = snakemake.params["data_configs"]
    cytotable_configs = snakemake.params["cytotable_config"]

    # executing conversion input file to parquet
    convert_to_parquet(
        input_file=plate_data,
        out_path=output_path,
        convert_configs=cytotable_configs,
        target_ext=general_configs,
    )


# executes the main function for conversion
if __name__ == "__main__":
    main()
