from pathlib import Path
from random import randint

if __name__ == "__main__":

    # create 5 random sqlite files
    for _ in range(5):
        rand_id = str(randint(1000, 9999))
        f_name = f"SQ{rand_id}.sqlite"

        with open(f_name, "w") as f:
            print(f"{f_name} created")

    # creating metadata dir
    path_obj = Path("metadata")
    path_obj.mkdir(exist_ok=True)

    # creating barcodes
    with open("barcodes.txt", "w") as f:
        print("bardcodes file created")

    # creating plate map file
    with open("plate_map.csv", "w") as f:
        print("plate map file created")
