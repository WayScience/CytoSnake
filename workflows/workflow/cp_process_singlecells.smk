"""
workflow: cp_process_singlecells.smk
"""


# importing modules
include: "../rules/common.smk"
include: "../rules/cytotable_convert.smk"
include: "../rules/normalize.smk"
include: "../rules/feature_select.smk"

# checking if configs if the users wants to generate consensus profiles
# if config["consensus"]:


#     include: "../rules/generate_consensus.smk"
rule all:
    input:
        CONVERTED_DATA_EXTENDED,
        NORMALIZED_DATA_EXPAND,
        SELECTED_FEATURE_DATA_EXPAND,
