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

```

##Usage
To run the phase diagram calculations, ensure that tc_python is installed and properly set up. Then, execute the script:

Clone or download the repository to your local machine.
Navigate to the directory containing the script.
Run the script using:
```bash
python phase_diagram_calculator.py
```
By default, the script uses the TCHEA6 Thermo-Calc database, computes phase diagrams over a temperature range of 500 K to 1200 K, and saves the results in a folder named results/. Each output file is named based on the composition it represents, such as Al0.05-Cr0.3-Co0.15-Fe0.2-Ni0.3.pkl.

To customize the calculation, modify the parameters in the script:

element_list: Define the elements in your system (e.g., ["Al", "Cr", "Co", "Fe", "Ni"]).
constant_C_of_others: Set fixed mole fractions for constant elements.
indices_of_changing_elements and indices_of_constant_elements: Specify which elements vary and which remain constant.
temperature_range: Define the temperature range for phase diagram calculations.
steps_in_temperature: Adjust the resolution of the temperature axis.
output_dir: Set the directory for saving results.
Results
The output files are stored as .pkl files, each corresponding to a specific composition. These files contain the phase diagram data grouped by stable phases. To load and analyze the data, use the included load_groups function:

```python
from phase_diagram_calculator import load_groups

data = load_groups("results/Al0.05-Cr0.3-Co0.15-Fe0.2-Ni0.3.pkl")
```
The data object contains structured information about the phase diagram, which can be processed or visualized. The script includes commented-out plotting code as a starting point for generating visualizations.

