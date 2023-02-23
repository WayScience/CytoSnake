from setuptools import find_packages, setup

setup(
    name="CytoSnake",
    version="0.0.1",
    url="https://github.com/WayScience/CytoPipe",
    author="Erik Serrano",
    packages=find_packages(),
    python_requires=">=3.10",
    entry_points={
        "console_scripts": [
            "cytosnake=cytosnake.cli.cmd:run_cmd",
        ]
    },
)
