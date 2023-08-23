# Functional Tests

- [Functional Tests](#functional-tests)
  - [Summary](#summary)
  - [Rational](#rational)
  - [Functional Testing Environment](#functional-testing-environment)
  - [Submitting Tests](#submitting-tests)
    - [Understanding the main steps of CytoSnake's Functional Tests](#understanding-the-main-steps-of-cytosnakes-functional-tests)
      - [Preparation](#preparation)
      - [Execution](#execution)
    - [Output Check](#output-check)
  - [Ending Remarks](#ending-remarks)

## Summary

Here, we provide you with `CytoSnake` functional test documentation.
It is designed to help you, as both a developer and user, gain a comprehensive understanding of the test cases and their expected outcomes.
By carefully reviewing this documentation, you can efficiently perform functional tests, validate `CytoSnake`'s behavior, and ensure its reliability in various scenarios.

## Rational

Functional tests are essential for evaluating `CytoSnake`'s behavior in accurately executing its intended functions with specific input parameters.
They primarily aim to validate the functionality of the program, ensuring it meets the specified requirements.
These tests involve feeding `CytoSnake`'s Command Line Interface (CLI) with various parameters to simulate user interactions and thoroughly monitor the user experience
By conducting functional tests, developers can ensure that all user-inputted parameters function correctly within `CytoSnake`, leading to a more reliable and user-friendly application.
Ultimately, the focus of these tests is to guarantee easy and simple user experience with `CytoSnake`.

## Functional Testing Environment

![functional_test_diagram](../../docs/images/functional-testing-diagram.svg)

> **Functional Testing Environment**: Schematic depicting `CytoSnake`'s functional testing environment.

Shown above is a basic illustration of the functional testing environment used by `CytoSnake`. 
This illustration introduces four fundamental components: the `CytoSnake` module(shown in blue-green), prepare_data function (shown in light-green), functional test module (shown in yellow), and fixture (shown in red) working together.

- **`CytoSnake module`**: is used to import `test_utils.py` module in order to use the `prepare_data()` function
- **`prepare_data()`**: allows developers to easily select data and transport them to the testing directory
- **`functional test`** (shown in yellow) contains the sets of instructions on how to test `CytoSnake`s'CLI. 
- **`fixture`**: Are [`pytest`](https://docs.pytest.org/en/6.2.x/fixture.html) fixtures that provide essential setup and teardown functionalities to make it easier for developers to form their tests.

## Submitting Tests

Adding functional tests is simple! 
You, as a developer, can create a new functional test within the `cytosnake/test/functional` directory and select a module where the test will be added.
This is similar to adding tests in traditional software packages, but using a consistent pattern to help maintainability.

Below, is an example to illustrate how to add a test in the `test_cli.py` module:

```python
# inside test_cli.py testing module
@pytest.mark.positive # or can be negative as well
def test_example(testing_dir) -> None:
    """
    [ Enter Documentation]
    test type: (positive or negative)
    rational: small description why this test is being done.
    input documentation: what inputs are being used in your test.
    """

    # step 1: prepare testing files and select input files
    datafiles = prepare_dataset(
        test_data_name="nf1-data", test_dir_path=testing_dir
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
To ensure that `pytest` recognizes your functions as tests, it's essential to begin the function name with `test_`.

To maintain clarity in your functional tests, please provide concise and well-structured documentation explaining the **test type** (positive or negative), **rationale**, and **inputs** (if any).
Here is a list below explaining the expected components of the documentation:

- **positive test**: indicates a successful execution with the provided input parameters.
- **negative test**: indicates that an error is supposed to be expected in the test.
- **rational**: reason behind this test, what is it trying to emulate?
- **inputs**: provide input documentation (e.g. explaining the dataset used).

Doing so will enable other developers in the `CytoSnake` community to fully comprehend the implemented functional tests.

### Understanding the main steps of CytoSnake's Functional Tests

The process of creating a functional test is straightforward, involving three key steps that developers must follow to build a successful and robust test.
These steps consist of:

- **preparation**: Preparing testing folder, dataset, and selecting inputs
- **execution**: Execute `CytoSnake` within testing older
- **assertion checks**: Checking desired outputs

#### Preparation

Preparing your dataset requires minimal effort with the help  `prepare_dataset` function.
This function readily prepares the dataset using just two parameters: the name of the dataset and the target directory for the test execution.

> **Note**: Thanks to the use of pytest's fixtures, the `testing_dir` variable already holds the path to the generated testing directory. As a result, users **are not required** to create their own testing directories when conducting functional tests.

The dataset names are derived from the folder names within the `cytosnake/tests/datasets` directory.
If the desired dataset is not available, developers can create their datasets; however, they must ensure they follow `CytoSnake's` input [requirements](https://cytosnake.readthedocs.io/en/latest/tutorial.html#modes).

After executing the `prepare_dataset()` function, the fixture generates an object called `DataFiles`, granting users full control over the selection of input files for `CytoSnake`.
For instance, if the chosen dataset comprises multiple plates, users can select a specific plate using Python indexing, as demonstrated in the example above.
Once the desired files are selected, users can proceed to prepare them for execution.

#### Execution

In the execution step, `CytoSnake` runs with the provided inputs within the testing directory.
Fortunately, the `prepare_dataset()` function not only generates the `DataFile` object but also informs the testing module about the execution location for `CytoSnake`.
This means after executing `prepare_dataset()`, the test automatically changes its directory to the testing folder, `testing_dir`, ensuring that all subsequent executions take place within this designated directory (look at figure).
This implicit change of directory simplifies the testing process and ensures that `CytoSnake` runs smoothly within the expected environment.

### Output Check

The output check step is for verifying if the expected outputs of your tests are generated.
As a user, you can perform checks for multiple outputs to ensure a passing test.
For instance, you can validate if the correct error is being raised by utilizing one of the test utility functions called `get_raised_error` (please refer to the documentation for more details).
This helps you to thoroughly examine different aspects of the test results and validate the accuracy of `CytoSnake`'s behavior.

## Ending Remarks

Congratulations! You now know to create functional tests effectively.
By understanding the rationale behind each step, following the straightforward process, and referring to the well-documented guidelines, you are fully equipped to create robust and reliable tests for `CytoSnake`.
So, go ahead and design some tests that not only ensure the stability of the user experience but also contribute to the growth of the `CytoSnake` community.

If you have any questions or need assistance with creating tests, don't hesitate to reach out to us.
Simply refer to the [issue](https://github.com/WayScience/CytoSnake/issues) form in the repository, and we'll be more than happy to provide help.
Your contributions are highly valued, and we want to ensure that your functional testing experience with `CytoSnake` is easy.
Feel free to ask any questions or add comments, and we'll be there to help you every step of the way!

Happy testing!
