from textwrap import dedent

class RunModeDoc:
    
    @property
    def __doc__(self):
        run_doc = 1

 
"""
Run Documentation

Cytopipe's run mode allows users to execute Cytopipes available workflows

USE CASE EXAMPLE:
-----------------
cytopipe run cp_process      # using default execute settings
cytopipe run cp_process -c 7 # executing workflow using at most 7 cores

Arguments
---------
help                Displays Run mode documentation 
workflow            Name of the workflow to execute
-c, --max_cores     Maximum number of cores to use for the workflow default=1
"""

"""_summary_
"""
        

