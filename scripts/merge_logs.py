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

    # convert list of logs into dataframe
    # sorting based on timestamp
    # convert sorted pandas df into list of log entries
    column_names = ["datetime", "debugg_type", "thread_id", "message"]
    log_df = pd.DataFrame(data=logs, columns=column_names)
    log_df["datetime"] = pd.to_datetime(log_df["datetime"], errors="coerce")
    log_df = log_df.sort_values("datetime", ascending=True)
    log_df["datetime"] = log_df["datetime"].astype(str)

    # write log file
    with open(outname, "w") as outlog:
        for row_data in log_df.values.tolist():
            log_str = " ".join(row_data) + "\n"
            outlog.write(log_str)

if __name__ in "__main__":

    # getting all logs and extracting all contents
    log_files = list(snakemake.input)
    outname = str(snakemake.output)

    logs = extract_logs(log_file_paths=log_files)
    write_merge_logs(logs=logs, outname=outname)




