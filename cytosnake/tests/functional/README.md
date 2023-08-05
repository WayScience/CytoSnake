# CytoSnake Functional Tests

Here, we provide you with `CytoSnake` functional test documentation.
It is designed to help you, as both a developer and user, to gain a comprehensive understanding of the test cases and their expected outcomes.
By carefully reviewing this documentation, you can efficiently perform functional tests, validate `CytoSnake`'s behavior, and ensure its reliability in various scenarios.

## Rational

Functional tests are essential for evaluating the performance of CytoSnake in accurately executing its intended functions with specific input parameters.
They primarily aim to validate the functionality of the program, ensuring it meets the specified requirements.
These tests involve feeding CytoSnake's Command Line Interface (CLI) with various parameters to simulate user interactions and thoroughly assess the user experience.
By conducting functional tests, developers can ensure that all user-inputted parameters function correctly within CytoSnake, leading to a more reliable and user-friendly application.
Ultimately, the focus of these tests is to guarantee a seamless and satisfactory user experience with CytoSnake.

## How to add tests

Adding functional tests is simple! You, as a developer, can create a new functional test within the `cytosnake/test/functional` directory and select a module where the test will be added.
Below, is an example to illustrate how to add a test in the `test_cli,py` module:

### Creating a tests

[ Add testing image here ]

```python
# inside test_cli.py testing module
def test_example(testing_dir) -> None:
    """
    [ Enter Documentation]
    test type: (positive or negative)
    rational:Small description why this test is being done.
    input documentation: what inputs are being used in your test.
    """

    # step1: prepare testing files and select input files
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

Above is an example of how a test should be written.
To ensure that `pytest` recognizes your functions as tests, it's essential to begin the function name with `test_.`"

To maintain clarity in your functional tests, please provide concise and well-structured documentation explaining the **test type** (positive or negative), **rationale**, and **inputs** (if any).
Doing so will enable other developers in the `CytoSnake` community to fully comprehend the implemented functional tests.

The process of creating a functional test is straightforward, involving four key steps that developers must follow to build a successful and robust test.
These steps consist of:

- **preparation**: Preparing testing folder and dataset
- **execution**: Execute `CytoSnake` within testing older
- **assertion checks**: Checking desired outputs

### Preparation

Preparing your dataset requires minimal effort with the help of the testing module's `prepare_dataset` function.
This function readily prepares the dataset using just two parameters: the name of the dataset and the target directory for the test execution.

> **Note**: Thanks to the use of pytest's fixtures, the testing_dir variable already holds the path to the generated testing directory. As a result, users are not required to creat their own testing directories when conducting functional tests.

The dataset names are derived from the folder names within the `cytosnake/tests/functional/dataset` directory.
If the desired dataset is not available, developers can create their own datasets; however, they must ensure they follow CytoSnake's input requirements.

After executing the `prepare_dataset` function, it generates an object called `DataFiles`, granting users full control over the selection of input files for `CytoSnake`.
For instance, if the chosen dataset comprises multiple plates, users can select a specific plate using Python indexing, as demonstrated in the example above.
Once the desired files are selected, users can proceed to prepare them for execution.

### Execution

In execution step, `CytoSnake` runs with the provided inputs within the testing directory.
Fortunately, the `prepare_dataset` function not only generates the `DataFile` object but also informs the testing module about the execution location for `CytoSnake`. This means after executing `prepare_dataset`, the test automatically changes its directory to the testing folder, `testing_dir`, ensuring that all subsequent executions take place within this designated directory. This implicit change of directory simplifies the testing process and ensures that CytoSnake runs smoothly within the expected environment.

### Assertion Check

The assertion check step is where if the expected outputs of you tests are generated. User can check for multiple outputs to ensure a passing test. For example, user can check if the correct error is being raised by using one of the test utils functions `get_raised_error` (refer to documentation )

The assertion check step is for verifying if the expected outputs of your tests are generated.
As a user, you can perform checks for multiple outputs to ensure a passing test.
For instance, you can validate if the correct error is being raised by utilizing one of the test utility functions called `get_raised_error` (please refer to the documentation for more details).
This helps you to thoroughly examine different aspects of the test results and validate the accuracy of CytoSnake's behavior.
