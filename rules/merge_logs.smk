# ----------------------------------------------------------------
# Documentation:
# Rule collects all generated logs from all porcessess and merges
# them into a single log file.
#
# individual log files are stored into an archive file along with
# the generated merged log.
#
# The archive file is taged with (Month-day-year)-(hour-min-sec)
# Example: 072922-083033_archived_logs
# ----------------------------------------------------------------

rule merge_logs:
    input:
        expand("logs/{logname}.log", logname=LOG_NAMES),
    output:
        "logs/CytoPipe_preprocess.log",
    script:
        "../scripts/merge_logs.py"
