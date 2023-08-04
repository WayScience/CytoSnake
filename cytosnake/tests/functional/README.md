# CytoSnake Functional Tests 

## Rational

Functional tests are essential for evaluating the performance of CytoSnake in accurately executing its intended functions with specific input parameters. 
They primarily aim to validate the functionality of the program, ensuring it meets the specified requirements. 
These tests involve feeding CytoSnake's Command Line Interface (CLI) with various parameters to simulate user interactions and thoroughly assess the user experience. 
By conducting functional tests, developers can ensure that all user-inputted parameters function correctly within CytoSnake, leading to a more reliable and user-friendly application. 
Ultimately, the focus of these tests is to guarantee a seamless and satisfactory user experience with CytoSnake.

## File structure


## How to add tests
Adding functional tests is simple! Developers are required to create a single function found within the `cytosnake/test/functional` directory and select a module where the tests is going to be added. 
Below is an example on how to add a tests.


### creating a tests 


```python
def test_example(testing_dir) -> None:
    """ 
    [ Enter Documentation]
    test type: (positive or negative)
    rational:Small description why this test is being done. 
    input documentation: what inputs are being used in your test.
    """

    # step1: prepare testing files
    datafiles = prepare_dataset(
        test_data_name="standard_sqlite_single", test_dir_path=testing_dir
    )

    # Selecting inputs
    plate = datafiles.plate_data[0]
    metadata = datafiles.metadata

    # Step 2: Execution of CytoSnake in testing directory
    cmd = f"cytosnake init -d {plate} -m {metadata}".split()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)

    # Step 3: Check if the test passed.
    assert proc.returncode == 0
```

Above is an example on how the functional tests appear. 
It takes in one parameter which is `testing_dir` a `pytest` fixture 

- talk about main function

Step 2 involves with execute the `CytoSnake` CLI command into the testing folder. 
Since the files of your selected testing data set has been trasnfered, 

- assertion


### add tests with dataset