import os
import glob
# assimptions:
# - must have a barcode map, metadata, and associated sql file
# - inputs must have consistent
#
#--------------------
# if multiple
# sql_files = golob.glob("input/*.sql")

# collecting all unqieuids from plate

rule aggregate:
    input:
        sql_files=expand("data/{plate_id}.sqlite", plate_id=PLATE_IDS),
        barcodes="data/barcode_platemap.csv",
        metadata="data/metadata"
    output:
        cell_counts=expand("results/preprocessing/{plate_id}.cell_counts.tsv", plate_id=PLATE_IDS),
        aggregate_profile=expand("results/preprocessing/{plate_id}.aggregate.csv.gz", plate_id=PLATE_IDS)
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/aggregate_cells.py"


rule annotate:
    input:
        barcodes="data/barcode_platemap.csv",
        aggregate_profiles=expand("results/preprocessing/{plate_id}.aggregate.csv.gz", plate_id=PLATE_IDS)
    output:
        expand("{plate_id}_augmented.csv.gz", plate_id=PLATE_IDS)
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/annotate.py"

























































# # NOTE: assuming one entry
# platename = "SQ00014613.sqlite".split(".")[0]
# rule aggergate_cells:
#     # NOTE: what are we assuming that user has?
#     # -- should we expect multiple?
#     input:
#         sql_file="SQ00014613.sqlite",
#         barcode="barcode_platemap.csv"

#     params:
#         strata="mage_Metadata_Plate Image_Metadata_Well"
#         imgage_cols="TableNumber ImageNumber Metadata_Site"
#         compression="gzip"
#     output:
#         count_out="results/{platename}.cell_counts.csv"
#         agger_out="data/{platename}.cell_aggregates.gz."
#         meta_out="metadata/{platename}.csv"
#     threads:
#     conda:
#         "../envs/cytominer_env.yaml"
#     # shell:
#         "echo scripts/aggergate.py -i {input.sql} -o {}