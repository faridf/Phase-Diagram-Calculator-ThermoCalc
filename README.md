# Phase-Diagram-Calculator-ThermoCalc
Phase Diagram Calculator in ThermoCalc with Python

This repository contains a Python script that automates the generation of phase diagram data using Thermo-Calc via the `tc_python` API. It systematically varies the compositions of a multi-component alloy system and computes the corresponding phase diagrams, saving the results for later analysis.

## Features

- Generate a mesh of compositions for a specified set of elements.
- Use Thermo-Calc to calculate phase diagrams at each composition point.
- Save phase diagram data as pickled Python dictionaries for easy retrieval and further processing.
- Customize the number of steps, temperature ranges, databases, and output directories.

## Requirements

- Python 3.7+ recommended.
- [Thermo-Calc Python API (tc_python)](https://thermocalc.com/products/thermo-calc/) and a valid license.
- NumPy
- Matplotlib (optional, if you plan to visualize results)
- pickle (standard library)

You can install Python dependencies (excluding Thermo-Calc which must be installed separately) using:
```bash
pip install -r requirements.txt
