import glob
from pathlib import Path
from datetime import datetime

import pandas as pd


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


def get_all_logs() -> list[str]:
    """Retuns all logs relative paths

    Returns
    -------
    list[str]
        relative paths of logs
    """
    log_rel_paths = glob.glob("./logs/*")
    log_abs_path = [str(Path(rel_path).absolute()) for rel_path in log_rel_paths]
    return log_abs_path


def write_merge_logs(logs: list[str], outname: str) -> None:
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
