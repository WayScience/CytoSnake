import shutil
from datetime import datetime
from pathlib import Path
from typing import Union

import pandas as pd


def archive_logs(files: Union[str, list[str]], dest: str) -> None:
    """Move single or list of files into

    Parameters
    ----------
    files : Union[str, list[str]]
        path of files to move
    dest : str
        Destinations to where the files will be move to
    """

    if isinstance(files, str):
        files = files.split()

    # checking if directory is passed
    if Path(dest).is_file():
        e_msg = "must be a path to a directory not to a file"
        raise TypeError(e_msg)
    elif not Path(dest).is_dir():
        print("Archive folder does not exists")

    # moving files
    for f_path in files:
        shutil.move(f_path, dest)


def extract_logs(log_file_paths: list[str]) -> list[str]:
    """Extracts all log contents within log files

    Parameters
    ----------
    log_file_paths : list[str]
        list of paths pointing to log files

    Return
    ------
    list[str]
        List of log contents
    """
    all_logs = []
    for log_file in log_file_paths:
        with open(log_file, "r") as logfile:
            for log in logfile:
                log = log.strip().replace("\n", "").split(" - ")
                all_logs.append(log)

    return all_logs


def sort_log_files(log_paths: list[str]) -> list[str]:
    """Sorts log files based on creation time

    Parameters
    ----------
    f_paths : list[str]
        list of all logs files

    Returns
    -------
    list[str]
        returns a list of log paths sorted based on creation
    """
    # collectime time created
    file_time = []
    for log_path in log_paths:
        log_path_obj = Path(log_path)
        mod_time = log_path_obj.stat().st_mtime
        time = datetime.fromtimestamp(mod_time)
        path_str = str(log_path_obj.absolute())
        file_time.append((time, path_str))

    file_time_df = pd.DataFrame(file_time, columns=["time_created", "path"])
    file_time_df["time_created"] = pd.to_datetime(file_time_df["time_created"])
    file_time_df = file_time_df.sort_values("time_created", ascending=True)
    sorted_paths = file_time_df["path"].tolist()
    return sorted_paths


def combine_logs(logs: list[str], outname: str) -> None:
    """Merges all logs into one file

    Parameters
    ----------
    logs : list[str]
        list of logged messages.

    outname : str
        Name of generated output file.

    Return
    ------
    None
        Creates merged log file
    """
    # sort logs by time created
    sorted_logs = sort_log_files(logs)

    # convert list of logs into dataframe
    # sorting based on timestamp
    # convert sorted pandas df into list of log entries
    column_names = ["datetime", "debugg_type", "thread_id", "message"]

    conts = extract_logs(sorted_logs)
    log_df = pd.DataFrame(data=conts, columns=column_names)

    # write log file
    with open(outname, "w") as outlog:
        for row_data in log_df.values.tolist():
            log_str = " ".join(row_data) + "\n"
            outlog.write(log_str)

    # archiving individual log files
    # -- using month-day-year and hour-min-sec as and id
    id_tag = datetime.now().strftime("%m%d%y-%H%M%S")
    dir_path = f"logs/{id_tag}-archived_logs"

    # -- creating directory for archiving logs
    archive_dir_obj = Path(dir_path)
    archive_dir_obj.mkdir(exist_ok=False)
    archive_dir_path = str(archive_dir_obj.absolute())

    archive_logs(files=sorted_logs, dest=archive_dir_path)
    shutil.copy(outname, archive_dir_path)


if __name__ == "__main__":

    # loading snakemake inputs
    log_files = list(snakemake.input)
    out_name = str(snakemake.output)

    # merge logs
    combine_logs(logs=log_files, outname=out_name)
