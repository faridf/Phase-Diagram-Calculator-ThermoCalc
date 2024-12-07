import numpy as np
import pickle
import os
from tc_python import *

def save_groups(groups, filename):
    """
    Save computed groups (phase diagram data) to a pickle file.

    Parameters
    ----------
    groups : dict
        Dictionary containing phase diagram data grouped by stable phases.
    filename : str
        The filename (including path) where data will be saved.
    """
    with open(filename, "wb") as f:
        pickle.dump(groups, f)

def load_groups(filename):
    """
    Load groups (phase diagram data) from a pickle file.

    Parameters
    ----------
    filename : str
        The filename (including path) from which to load the data.

    Returns
    -------
    dict
        Dictionary containing phase diagram data grouped by stable phases.
    """
    with open(filename, "rb") as f:
        groups = pickle.load(f)
    return groups

def mesh_concentrations(element_list, n, indices_of_changing_elements,
                        indices_of_constant_elements, constant_concentration_of_other_elements):
    """
    Generate a matrix of input concentrations for a system of elements.

    The function creates a series of compositions by varying two elements 
    (as specified by indices_of_changing_elements) while keeping the other 
    elements (indices_of_constant_elements) constant.

    Parameters
    ----------
    element_list : list of str
        List of element symbols, e.g. ["Al", "Cr", "Co", "Fe", "Ni"].
    n : int
        Number of iterations for generating the concentration mesh.
    indices_of_changing_elements : list of int
        Indices of elements in element_list that will be varied.
    indices_of_constant_elements : list of int
        Indices of elements in element_list that remain constant.
    constant_concentration_of_other_elements : float
        The fixed concentration assigned to each of the constant elements.

    Returns
    -------
    np.ndarray
        A 2D array of shape (len(element_list), n) containing the mole 
        fractions for each element at each iteration.
    """
    n_elements = len(element_list)
    input_concentrations = np.zeros((n_elements, n))
    mole_fractions = np.linspace(0, 1 - 3 * constant_concentration_of_other_elements, n)
    
    for i in range(n):
        # Set constant elements
        for idx in indices_of_constant_elements:
            input_concentrations[idx, i] = constant_concentration_of_other_elements
        # Varying elements
        # First changing element
        input_concentrations[indices_of_changing_elements[0], i] = mole_fractions[i]
        # Second changing element
        input_concentrations[indices_of_changing_elements[1], i] = (
            1 - mole_fractions[i] - 3 * constant_concentration_of_other_elements
        )
    return input_concentrations


def run_phase_diagram_calculations(element_list, constant_C_of_others, n=15, 
                                   indices_of_changing_elements=[0, 1], 
                                   indices_of_constant_elements=[2, 3, 4],
                                   database="TCHEA6", 
                                   temperature_range=(500,1200), 
                                   steps_in_temperature=60,
                                   output_dir="results"):
    """
    Run a series of Thermo-Calc phase diagram calculations over a grid of compositions.

    Parameters
    ----------
    element_list : list of str
        List of element symbols, e.g. ["Al", "Cr", "Co", "Fe", "Ni"].
    constant_C_of_others : list of float
        List of constant mole fractions to assign to the "other" (constant) elements.
    n : int, optional
        Number of iterations for varying concentrations along the composition axis, by default 15.
    indices_of_changing_elements : list of int, optional
        Indices of elements in element_list that will be varied, by default [0, 1].
    indices_of_constant_elements : list of int, optional
        Indices of elements in element_list that remain constant, by default [2, 3, 4].
    database : str, optional
        Name of the Thermo-Calc database to use, by default "TCHEA6".
    temperature_range : tuple, optional
        Min and max temperature (K), by default (500, 1200).
    steps_in_temperature : int, optional
        Number of steps in temperature to discretize the calculation, by default 60.
    output_dir : str, optional
        Directory to store the results, by default "results".

    Notes
    -----
    This function requires Thermo-Calc Python API (tc_python) and a valid license.
    The `TCPython()` context manager is used to handle sessions.

    Results are saved as pickle files named using the composition 
    of the system at each calculation point.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    counter = 0
    for j in constant_C_of_others:
        input_concentrations = mesh_concentrations(
            element_list=element_list,
            n=n, 
            indices_of_changing_elements=indices_of_changing_elements, 
            indices_of_constant_elements=indices_of_constant_elements, 
            constant_concentration_of_other_elements=float(j)
        )

        # We skip first and last indices to avoid trivial endpoints as in original code
        for i in range(1, input_concentrations.shape[1] - 1):
            counter += 1
            print(f"Calculating system #{counter}")

            # Extract mole fractions for current system
            Al_x = float(input_concentrations[0, i])
            Cr_x = float(input_concentrations[1, i])
            Co_x = float(input_concentrations[2, i])
            Fe_x = float(input_concentrations[3, i])
            Ni_x = float(input_concentrations[4, i])

            # Determine maximum concentration for the first axis (Al in this case)
            max_concentration_on_x_axis = 1.0 - Cr_x - Co_x - Fe_x - Ni_x
            max_concentration_on_x_axis = round(max_concentration_on_x_axis, 3)

            # Run the Thermo-Calc calculation
            try:
                with TCPython() as session:
                    calculator = (
                        session
                        .select_database_and_elements(database, element_list)
                        .get_system()
                        .with_phase_diagram_calculation()
                        .with_first_axis(
                            CalculationAxis(ThermodynamicQuantity.mole_fraction_of_a_component("Cr"))
                            .set_min(0)
                            .set_max(max_concentration_on_x_axis)
                            .with_axis_type(Linear().set_max_step_size(0.025 / max_concentration_on_x_axis))
                        )
                        .with_second_axis(
                            CalculationAxis(ThermodynamicQuantity.temperature())
                            .set_min(temperature_range[0])
                            .set_max(temperature_range[1])
                            .with_axis_type(AxisType.linear().set_min_nr_of_steps(steps_in_temperature))
                        )
                        .enable_global_minimization()
                        .set_condition(ThermodynamicQuantity.temperature(), 1000)
                        .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("Cr"), Cr_x)
                        .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("Co"), Co_x)
                        .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("Fe"), Fe_x)
                        .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("Ni"), Ni_x)
                    )

                    phase_diagram = calculator.calculate(timeout_in_minutes=15)
                    phase_diagram = (
                        phase_diagram
                        .add_coordinate_for_phase_label(0.5, 2000)
                        .get_values_grouped_by_stable_phases_of(
                            ThermodynamicQuantity.mole_fraction_of_a_component("Al"),
                            ThermodynamicQuantity.temperature()
                        )
                    )

                    # Save the results
                    filename = f"Al{Al_x}-Cr{Cr_x}-Co{Co_x}-Fe{Fe_x}-Ni{Ni_x}.pkl"
                    filepath = os.path.join(output_dir, filename)
                    save_groups(phase_diagram, filepath)
                    print(f"Saved: {filepath}")

            except UnrecoverableCalculationException:
                print("Could not calculate. Continuing with next...")

if __name__ == "__main__":
    # Example usage:
    element_list = ["Al", "Cr", "Co", "Fe", "Ni"]
    constant_C_of_others = [0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3]
    run_phase_diagram_calculations(element_list, constant_C_of_others)
