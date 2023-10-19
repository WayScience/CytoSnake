# Benchmarking Workflows

## Introduction to Benchmarking with Cytosnake

Cytosnake is a powerful Python package known for its extensive collection of workflows designed for analyzing image-based profiles. To ensure these workflows perform at their best, it's essential to conduct benchmarking. Benchmarking allows us to assess and optimize the efficiency and reliability of Cytosnake's image analysis capabilities. This tutorial will guide you through the benchmarking process, helping you gauge and improve the performance of these image-based profiling workflows. Whether you're a developer seeking to enhance Cytosnake's capabilities or a user interested in its performance, this guide will equip you with the necessary knowledge.

Here in Cytosnake, we use Memray to benchmark our workflows, making it a great tool for tackling memory usage issues, identifying memory leaks, and pinpointing code hotspots that result in excessive allocations.
Memray's notable features include tracing every function call, handling native calls in C/C++ libraries, minimal application slowdown during profiling, generating diverse reports like flame graphs, compatibility with Python threads, and support for native threads.

In this guide, you'll gain an understanding of how to conduct benchmarks with Cytosnake.
You'll learn not only on how to execute benchmarking but also how to interpret and make the most of the benchmarking outputs.

## Enable benchmarking in cytosnake

To enable benchmarking, you can configure it by navigating to the `config/configurational.yaml` file and opening it using your preferred text editor.
 Within the file, set the `enable_profiling` to `True` like this:

```yaml
enable_profiling: True
```

By doing so, you inform CytoSnake that you intend to perform benchmarking on the workflow you are about to utilize.

## Executing benchmarking

This steps assumes that you've already initialized your data using the `cytosnake init` mode. To execute benchmarking, the process is as the same as running a typical workflow:

```bash
cytosnake run cp_process
```

Behind the scenes, benchmarking will be automatically initialized, and it will start benchmarking each step within the workflow.
