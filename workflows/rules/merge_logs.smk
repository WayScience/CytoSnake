"""
rule module: merge_logs.smk

Collects all log files generated within each rule module and merges it into
one log file

The log file is tagged with (Month-day-year)-(hour-min-sec)
Example: 072922-083033_archived_logs

Paramters:
Inputs:
    No user defined outputs, searches individual logs in the `logs/` folder
Output:
    Merged log file


Returns
    Merged log file stored in the `logs/` directory
"""


rule merge_logs:
    input:
        expand("logs/{logname}.log", logname=LOG_NAMES),
    output:
        "logs/CytoSnake_preprocess.log",
    script:
        "../scripts/merge_logs.py"
