import metaspread.configs
import pandas as pd
import shutil
import mesa
import ast
import os
from pathlib import Path
import random
from datetime import datetime

# To run this code you must be in the parent folder of the program

def save_configs(simulations_dir, new_simulation_folder, config_var_names, max_steps, data_collection_period):
    # Saves the simulation configuration
    print(f"\t Saving all the simulations parameters at: {os.path.join(simulations_dir, new_simulation_folder, 'configs.csv')}")
    values = [getattr(metaspread.configs, i) for i in config_var_names]
    names = config_var_names

    #add configurations that are not in the global variables
    names += ['max_steps', 'data_collection_period']
    values += [max_steps, data_collection_period]
    df_vars = pd.DataFrame({"Names": names, "Values": values})
    df_vars = df_vars.set_index("Names")
    path = os.path.join(simulations_dir, new_simulation_folder, 'configs.csv')
    df_vars.to_csv(path)


def run_simulation(simulation_id, max_steps, data_collection_period, save_path=Path("."), loaded_simulation_path=""):
    # n = random.randint(1, 100)
    # load configs file from a previous simulation or loads the general configs file
    # print(loaded_simulation_path)
    loaded_simulation_path= loaded_simulation_path.strip('\"')
    if loaded_simulation_path != "":
        configs_path = os.path.join(loaded_simulation_path, "configs.csv")
        config_var_names = metaspread.configs.load_simulation_configs_for_reloaded_simulation(configs_path)
    else:
        package_dir = os.path.dirname(metaspread.__file__)
        configs_path = os.path.join(package_dir, "simulations_configs.csv")
        # print(f"simrunner.py configs_path: {configs_path}")
        config_var_names = metaspread.configs.init_simulation_configs(configs_path)
    
    # Parameters for this simulation
    number_of_initial_cells = metaspread.configs.number_of_initial_cells # Number of cancer cells
    gridsize     = metaspread.configs.gridsize
    grids_number = metaspread.configs.grids_number
    width        = gridsize
    height       = gridsize

    # Name of the directories
    simulations_dir = save_path / "Simulations"
    os.makedirs(simulations_dir, exist_ok=True)
    if loaded_simulation_path != "":
        cells_path = os.path.join(loaded_simulation_path, "CellsData.csv")
        df = pd.read_csv(cells_path)
        loaded_max_step = max(df["Step"])
        new_simulation_folder = os.path.normpath(loaded_simulation_path)
        new_simulation_folder = os.path.basename(new_simulation_folder)
        new_simulation_path = os.path.join(simulations_dir, new_simulation_folder)
    else:
        df = pd.DataFrame()
        loaded_max_step = 0
        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        new_simulation_folder = f"Sim-{simulation_id}-Date-{current_time}"

        # Creates the path for the new simulation
        new_simulation_path = os.path.join(simulations_dir, new_simulation_folder)
        pathMmp2 = os.path.join(new_simulation_path, "Mmp2")
        pathEcm = os.path.join(new_simulation_path, "Ecm")
        pathVasculature = os.path.join(new_simulation_path, "Vasculature")
        pathTimeOfPopulation = os.path.join(new_simulation_path, "Time when grids were populated")

        # Create folder for all cells analysis, for Mmp2 matrices and Ecm matrices
        if not os.path.exists(new_simulation_path):
            print(f'\t Folder for this simulation: {new_simulation_path}')
            print(f'\t Saving agents data at: {new_simulation_path}')
            print(f'\t Saving Mmp2 data at: {pathMmp2}')
            print(f'\t Saving Ecm data at: {pathEcm}')
            print(f'\t Saving Vasculature data at: {pathVasculature}')

            os.makedirs(new_simulation_path)
            os.makedirs(pathMmp2)
            os.makedirs(pathEcm)
            os.makedirs(pathVasculature)
            os.makedirs(pathTimeOfPopulation)
        # If there is already a simulation you skip it
        else:
            return print("This simulation already exists!")

    # Run the simulation and saves the data
    save_configs(simulations_dir, new_simulation_folder, config_var_names, max_steps, data_collection_period)
    model = metaspread.CancerModel(
        number_of_initial_cells,
        width,
        height,
        grids_number,
        max_steps,
        data_collection_period,
        new_simulation_path,
        loaded_simulation_path)
    for i in range(max_steps):
        model.step()
    print(f'Finished the simulation at time step {model.schedule.time}!')
    return model