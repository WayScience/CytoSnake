from setuptools import setup, find_packages

setup(
    name="cytopipe",
    version="0.0.1",
    url="https://github.com/WayScience/CytoPipe",
    author="Way Lab",
    packages=find_packages(),
    python_requires=">=3.10",
    entry_points={
        "console_scripts": [
            "cytopipe=cytopipe.cli.cmd:run_cmd",
        ]
    },
)
