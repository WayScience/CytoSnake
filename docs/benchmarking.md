# Benchmarking Workflows

<!-- used for displaying doctree when selecting file -->
```{toctree}
:maxdepth: 4
```

## Introduction to Benchmarking with Cytosnake

Cytosnake is a Python package known for its collection of workflows designed for performing image-based profiling.
To ensure these workflows perform at their best, it's essential to conduct benchmarking.
Benchmarking allows us to assess and optimize the efficiency and reliability of Cytosnake's image analysis processes.
This tutorial will guide you through the benchmarking process, helping you gauge and improve the performance of these image-based profiling workflows.
Whether you're a developer seeking to enhance Cytosnake's functionalities or a user interested in its performance, this guide will equip you with the necessary knowledge.

Here in Cytosnake, we use Memray to benchmark our workflows, making it a great tool for tackling memory usage issues, identifying memory leaks, and pinpointing code hotspots that result in excessive allocations.
Memray's notable features include tracing every function call, handling native calls in C/C++ libraries, minimal application slowdown during profiling, diverse report generation like flame graphs, compatibility with Python threads, and support for native threads.

In this guide, you'll gain an understanding of how to conduct benchmarks with Cytosnake.
You'll learn not only on how to execute benchmarking but also how to interpret and make the most of the benchmarking outputs.

## Enable benchmarking in cytosnake

To enable benchmarking, modify the `config/configurational.yaml` file and opening it using your preferred text editor.
Once inside the config file, set the `enable_memory_tracking` to `True`, informing CytoSnake that you intend to perform benchmarking on the workflow.

```yaml
enable_profiling: True
```

## Executing benchmarking

This steps assumes that you've already initialized your data using the `cytosnake init` mode.
To execute benchmarking, the process is as the same as running a typical workflow:

```bash
cytosnake run cp_process
```

Behind the scenes, benchmarking will be automatically initialized, and it will start benchmarking each step within the workflow.


## Examining the Benchmarks

Upon completion of the benchmarks, a `benchmarks/` folder will be generated in your `ProjectDirectory`. To understand more about the `ProjectDirectory`, refer to CytoSnake's documentation available [here](https://cytosnake.readthedocs.io/en/latest/tutorial.html#setting-up-files).

Two methods are available for analyzing the generated benchmark files:

### Method 1: Using Memray

You can employ the `memray` command-line tool to extract information from the binarized files and convert them into `.json` files. Execute the following command:

```bash
memray stats <BINFILE> --json --output <JSON_OUT_NAME>
```

For detailed instructions on extracting binarized `memray` benchmark outputs, visit [this link](https://bloomberg.github.io/memray/stats.html). The resulting JSON files contain benchmarking information captured within the workflow.

### Method 2: Using CytoSnake-Benchmark

Another approach to analyze the generated `benchmarks/` folder is by utilizing the dedicated benchmarking directory available at [CytoSnake-Benchmarks](https://github.com/WayScience/CytoSnake-Benchmarks).

Requirements include the `benchmarks/` folder and the path where all your data is stored. If you have executed `cytosnake run`, then the path you provide should be the `data` folder. Refer to the documentation there for guidance on conducting benchmarking.