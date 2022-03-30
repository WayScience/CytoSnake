
def test2():
    """Tests if snakemake variables are global or not
    """

    print("Testing function output")
    output = dict(snakemake.output)
    print(output)

if __name__ in "__main__":

    print("reading inputs")

    print(snakemake.input.sql_file)
    print(snakemake.input.barcode)
    print(snakemake.input.metadata)

    output = dict(snakemake.output)
    with open(str(output["cell_counts"]), "w") as f:
            f.write("output1")


    with open(str(output["aggregate_profile"]), "w") as f:
            f.write("output1")


    print("Testing is snakemake variables are global or not")
    test2() # yes it s! 

    print("testing daatypr")
    print(type(snakemake.output["cell_counts"]))