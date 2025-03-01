import mesa
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import os
import json
import ast
from metaspread.cancercell import CancerCell
from metaspread.vessel import Vessel
from metaspread.quasicircle import find_quasi_circle
from matplotlib import pyplot as plt
from matplotlib import cm
# from Classes.configs import *
import metaspread.configs
# import pickle


def get_cluster_survival_probability(cluster):
    """
    Takes in a tuple representing a cluster, returns the survival probabiltiy.
    
    Input:
        Cluster: a tuple representing the cancer cells cluster, where the first 
        element corresponds to the amount of Mesenchymal cells, and the second
        corresponds to the amount of Epithelial cells.
    Returns:
        The corresponding survival probability, according to if its a songle-cell
        cluster or a multi-cellular one.
    """
    if cluster[0] < 0:
        raise Exception(f"Error! Mesenchymal cells are negative: {cluster[0]}")
    if cluster[1] < 0:
        raise Exception(f"Error! Epithelial cells are negative: {cluster[1]}")
    if sum(cluster) == 1:
        return (single_cell_survival)
    elif sum(cluster) > 1:
        return (cluster_survival)
    elif sum(cluster) == 0:
        raise Exception(f"Error, no cells in cluster!")
    else:
        raise Exception(f"Error, nothing returned for cluster survival probability" )
    

def count_total_cells(model):
    """Counts all the cells present in the model, in all sites.

    Input:
        model: CancerModel object.
    Returns:
        amount_of_cells (int): the total amount of cells in every site,
        NOT considering the vasculature    
    """
    return sum(model.cancer_cells_counter)

def count_vasculature_cells(model):
    """"
    Counts the total amount of cells in the vasculature

    Input: 
        model: CancerModel object.
    Returns:
        amount_of_cells (int): the total amount of cells in the vasculature
    """
    amount_of_cells = sum([sum(x+y for x,y in value) for value in model.vasculature.values()])
    return amount_of_cells

class CancerModel(mesa.Model):
    """
    Class for the model.

    Attributes:
    ---------------
    number_of_initial_cells: int
        number of initial cancer cells
    width: int
        the width of each grid
    height: int
        the height of each grid
    grids_number: int
        the amount of sites, consideiring the initial tumour site + the secondary sites
    seed: int
        the seed used for the random number generation for all the simulation.
        If None, the random one will be selected by default.

    Methods:
    ---------------
    _initialize_grids__()
        initialize the grid with the initial bessel and cancer cell population
    proliferate(cell_type)
        Duplicates every cancer cell in the model of the cell_type phenotype
    calculate_environments(mmp2, ecm)
        Calculates the next step for the given arrays of mmp2 and ecm concentrations
    disaggregate_clusters(time)
        For a given time, it will dissagregate single cells from clusters
    """

    def __init__(self, number_of_initial_cells, width, height, grids_number, max_steps, data_collection_period, new_simulation_folder, loaded_simulation_path="", fixed_p_left=None, fixed_p_right=None, fixed_p_top=None, fixed_p_bottom=None, seed=None):
        super().__init__()  
        # self.simulations_dir = "Simulations"
        
        self.fixed_p_left=fixed_p_left
        self.fixed_p_right=fixed_p_right
        self.fixed_p_top=fixed_p_top
        self.fixed_p_bottom=fixed_p_bottom
        self.vasculature = {}
        self.number_of_initial_cells = number_of_initial_cells
        self.width = width
        self.height = height
        self.phenotypes = ["mesenchymal", "epithelial"]
        self.grid_vessels_positions = [[] for _ in range(grids_number)]
        self.current_agent_id = 0
        self.max_steps = max_steps
        self.data_collection_period = data_collection_period
        self.new_simulation_folder  = new_simulation_folder
        self.mesenchymal_count = [np.zeros((width, height), dtype=float) for _ in range(grids_number)]
        self.epithelial_count = [np.zeros((width, height), dtype=float) for _ in range(grids_number)]
        self.grids_number = grids_number
        self.grids = [mesa.space.MultiGrid(width, height, False) for _ in range(self.grids_number)]
        self.grid_ids = [i+1 for i in range(self.grids_number)]
        self.cancer_cells_counter = [0] * grids_number
        self.time_grid_got_populated = [-1 for _ in range(self.grids_number)]
        self.schedule = mesa.time.RandomActivation(self)
        #list of numpy arrays, representing mmp2 and ecm concentration in each grid
        self.mmp2 = [np.zeros((2, width, height), dtype=float) for _ in range(grids_number)]
        self.ecm = [np.ones((2, width, height), dtype=float) for _ in range(grids_number)]
        self.loaded_max_step = 0
        self.previous_cell_data = pd.DataFrame()

        if loaded_simulation_path != "":
            print(f"Loading simulation at {loaded_simulation_path}!")
            configs_path = os.path.join(loaded_simulation_path, "configs.csv")
            config_var_names = metaspread.configs.load_simulation_configs_for_reloaded_simulation(configs_path)
            #load the configs so we can use them in this module as globals
            for var in config_var_names:
                globals()[var] = getattr(metaspread.configs, var)
            self.load_previous_simulation(loaded_simulation_path)
            self.previous_cell_data = pd.read_csv(os.path.join(loaded_simulation_path, "CellsData.csv"), index_col=0)
        else:
            print("Starting simulation from zero!")
            configs_path = "simulations_configs.csv"
            config_var_names = metaspread.configs.init_simulation_configs(configs_path)
            #load the configs so we can use them in this module as globals
            for var in config_var_names:
                globals()[var] = getattr(metaspread.configs, var)
            self._initialize_grids()
            self.doubling_time_counter_M = doubling_time_M
            self.doubling_time_counter_E = doubling_time_E
        self.datacollector = mesa.DataCollector(
            model_reporters={"Total cells": count_total_cells}, agent_reporters={"Position": "pos", "Agent Type": "agent_type", "Phenotype": "phenotype", "Ruptured": "ruptured", "Grid": "grid_id"})

    def step(self):
        """Advance the model by one step.
        
        The step function of the model will be called on each of the siulations steps.
        It will correctly allocate the cells coming from the vaculature, if any.
        Then, it will calculate the ECM and MMP2 concentration changes for this step,
        proliferate the cells if its due, and collect the data if it is the appropiate time.

        Input: none
        Returns: none
        """       
        if self.schedule.time in self.vasculature: # Add keys
            self.disaggregate_clusters(self.schedule.time)
            surviving_clusters = [cluster for cluster in self.vasculature[self.schedule.time] if self.random.random() < get_cluster_survival_probability(cluster)]
            del self.vasculature[self.schedule.time]
            for cluster in surviving_clusters:
                selected_site = self.random.choices(range(1,self.grids_number), weights=extravasation_probs[0:self.grids_number-1])[0]
                arriving_point = self.random.choice(self.grid_vessels_positions[selected_site])
                x,y = arriving_point
                on_left_border    = self.grids[selected_site].out_of_bounds((x-1,y))
                on_right_border   = self.grids[selected_site].out_of_bounds((x+1,y))
                on_top_border     = self.grids[selected_site].out_of_bounds((x,y+1))
                on_bottom_border  = self.grids[selected_site].out_of_bounds((x,y-1))
                possible_places = self.grids[selected_site].get_neighborhood(arriving_point, moore=False, include_center=False)
                number_of_ccells_in_arriving_point ={}
                for x2,y2 in possible_places:
                    number_of_ccells_in_arriving_point[x2,y2] = len([agent for agent in self.grids[selected_site].get_cell_list_contents([(x2,y2)]) if agent.agent_type == "cell"])
                for tuple_index, ccells_amount in enumerate(cluster):
                    cell_type = "mesenchymal" if tuple_index == 0 else "epithelial"
                    while ccells_amount > 0:
                        if not on_left_border and carrying_capacity > number_of_ccells_in_arriving_point[x-1,y]:
                            ccell = CancerCell(self.current_agent_id, self, self.grids[selected_site], self.grid_ids[selected_site], cell_type, self.ecm[selected_site], self.mmp2[selected_site])
                            self.current_agent_id += 1
                            self.grids[selected_site].place_agent(ccell, (x-1,y)) 
                            number_of_ccells_in_arriving_point[x-1,y] += 1
                            self.cancer_cells_counter[selected_site] += 1
                            self.schedule.add(ccell)
                        elif not on_right_border and carrying_capacity > number_of_ccells_in_arriving_point[x+1,y]:
                            ccell = CancerCell(self.current_agent_id, self, self.grids[selected_site], self.grid_ids[selected_site], cell_type, self.ecm[selected_site], self.mmp2[selected_site])
                            self.current_agent_id += 1
                            self.grids[selected_site].place_agent(ccell, (x+1,y))
                            number_of_ccells_in_arriving_point[x+1,y] += 1
                            self.cancer_cells_counter[selected_site] += 1
                            self.schedule.add(ccell)
                        elif not on_bottom_border and carrying_capacity > number_of_ccells_in_arriving_point[x,y-1]:
                            ccell = CancerCell(self.current_agent_id, self, self.grids[selected_site], self.grid_ids[selected_site], cell_type, self.ecm[selected_site], self.mmp2[selected_site])
                            self.current_agent_id += 1
                            self.grids[selected_site].place_agent(ccell, (x,y-1))
                            number_of_ccells_in_arriving_point[x,y-1] += 1
                            self.cancer_cells_counter[selected_site] += 1
                            self.schedule.add(ccell)
                        elif not on_top_border and carrying_capacity > number_of_ccells_in_arriving_point[x,y+1]:
                            ccell = CancerCell(self.current_agent_id, self, self.grids[selected_site], self.grid_ids[selected_site], cell_type, self.ecm[selected_site], self.mmp2[selected_site])
                            self.current_agent_id += 1
                            self.grids[selected_site].place_agent(ccell, (x,y+1))
                            number_of_ccells_in_arriving_point[x,y+1] += 1
                            self.cancer_cells_counter[selected_site] += 1
                            self.schedule.add(ccell)
                        ccells_amount -= 1
                    
        #Perform ECM and MMP2 calculations
        self.calculate_environment(self.mmp2, self.ecm)
        
        # Proliferation
        # Counters are used so when loading a simulation the behaviour does not change, compared to use self.schedule.time % doubling_time_M == 0
        if (self.doubling_time_counter_M == 0 and self.schedule.time != 0):
            self.proliferate("mesenchymal")
            self.doubling_time_counter_M = doubling_time_M

        if (self.doubling_time_counter_E == 0 and self.schedule.time != 0):
            self.proliferate("epithelial")
            self.doubling_time_counter_E = doubling_time_E
                
        self.doubling_time_counter_E -= 1
        self.doubling_time_counter_M -= 1

        print(f'Step number: {self.schedule.time + self.loaded_max_step}', end="")
        print("\r", end="")
        self.schedule.step()
        
        #At the end of each step, check if the grid has been populated, and if it happened, store the time step when it did
        for index, time in enumerate(self.time_grid_got_populated):
            if time == -1: #if it has not been populated already, we check:
                if self.cancer_cells_counter[index] > 0:
                    self.time_grid_got_populated[index] = self.schedule.time + self.loaded_max_step

        # Saving of non agents data every data_collection_interval steps
        if (self.schedule.time != 0 and (self.schedule.time % self.data_collection_period == 0)) \
            or self.schedule.time == self.max_steps:
            self.datacollector.collect(self)
            current_agents_state = self.datacollector.get_agent_vars_dataframe()
            current_agents_state = current_agents_state.reset_index(level=["Step", "AgentID"])
            path_to_save = os.path.join(self.new_simulation_folder, f'CellsData.csv')
            if not self.previous_cell_data.empty:
                current_agents_state["Step"] += self.loaded_max_step
                current_agents_state = pd.concat([self.previous_cell_data, current_agents_state])
            current_agents_state.to_csv(path_to_save)
            #pickling a model could be an option in the future
            # backup_file_path = os.path.join(self.new_simulation_folder, "Backup", "backup.p")
            # with open(backup_file_path, "wb") as f:
            #     pickle.dump(self, f)
            df_time_grids_got_populated = pd.DataFrame()
            for grid_id in self.grid_ids:
                new_mmp2_df = pd.DataFrame(self.mmp2[grid_id-1][0,:,:])
                mmp2CsvName = f"Mmp2-{grid_id}grid-{self.schedule.time + self.loaded_max_step}step.csv"
                path_to_save = os.path.join(self.new_simulation_folder, "Mmp2", mmp2CsvName)
                new_mmp2_df.to_csv(path_to_save)

                new_ecm_df = pd.DataFrame(self.ecm[grid_id-1][0,:,:])
                EcmCsvName = f"Ecm-{grid_id}grid-{self.schedule.time + self.loaded_max_step}step.csv"
                path_to_save = os.path.join(self.new_simulation_folder, "Ecm", EcmCsvName)
                new_ecm_df.to_csv(path_to_save)

                df_time_grids_got_populated[f"Time when grid {grid_id} was first populated"] = [self.time_grid_got_populated[grid_id-1]]
                df_time_grids_got_populated_csv_name = f"Cells-are-present-grid-{grid_id}-{self.schedule.time + self.loaded_max_step}step.csv"
            path_to_save = os.path.join(self.new_simulation_folder, "Time when grids were populated", df_time_grids_got_populated_csv_name)
            df_time_grids_got_populated.to_csv(path_to_save)

            # Saves vasculature data
            # {key: list of clusters} -> {timestep: [(number of Mcells, number of Ecells), ..., (..., ...)]}
            vasculature_json = json.dumps(self.vasculature)
            
            vasculature_json_name = f"Vasculature-{self.schedule.time + self.loaded_max_step}step.json"
            path_to_save = os.path.join(self.new_simulation_folder, "Vasculature", vasculature_json_name)
            
            with open(path_to_save, 'w') as f:
                f.write(vasculature_json)
                
            # Saves cancer cells data as a backup in case the simulation fails
            # _, current_model_data = mesa.batchrunner._collect_data(self, self.data_collection_period-1)
            # df_current_model_data = pd.DataFrame(current_model_data)
            # df_current_model_data["Step"] = self.data_collection_period
            # for step in range(self.data_collection_period * 2, self.schedule.time, self.data_collection_period):
            #     _, step_model_data = mesa.batchrunner._collect_data(self, step-1)
            #     df_step_model_data = pd.DataFrame(step_model_data)
            #     df_step_model_data["Step"] = step + self.loaded_max_step
            #     df_current_model_data = pd.concat([df_current_model_data, df_step_model_data])
            # path_to_save = os.path.join(self.new_simulation_folder, f'CellsData.csv')
            # df_current_model_data.to_csv(path_to_save)



    def proliferate(self, cell_type):
        """"
        Duplicates every cell of cell_type phenotype in every site of the model

        Input: none
        Returns: none
        """
        for agent in self.schedule.agents:
            if agent.agent_type == "cell":
                x, y = agent.pos
                amount_of_cells = len([cell for cell in agent.grid.get_cell_list_contents([(x, y)]) if cell.agent_type == "cell"])
                if carrying_capacity > amount_of_cells and agent.phenotype == cell_type:
                    # print("Created new cell!!")
                    new_cell = CancerCell(self.current_agent_id, self, agent.grid, agent.grid_id, agent.phenotype, agent.ecm, agent.mmp2)
                    self.current_agent_id += 1
                    self.schedule.add(new_cell)
                    agent.grid.place_agent(new_cell, (x,y))
                    self.cancer_cells_counter[agent.grid_id - 1] += 1
        


    def load_previous_simulation(self, path_to_simulation):
        """
        Loads the last step of a previously computed simulation as the initial condition of this model

        Input: The simulation's path to be loaded
        Returns: none
        """
        #load mmp2 and ecm 
        mmp2_files_path = os.path.join(path_to_simulation, "Mmp2")
        ecm_files_path  = os.path.join(path_to_simulation, "Ecm")
        for grid_number in range(self.grids_number):
            mmp2_files = [file for file in os.listdir(mmp2_files_path) if file.startswith(f"Mmp2-{grid_number+1}grid")]
            ecm_files  = [file for file in os.listdir(ecm_files_path) if file.startswith(f"Ecm-{grid_number+1}grid")]
            mmp2_files.sort(key = lambda file_name: int(file_name.split('step')[0][11:]))
            ecm_files.sort(key  = lambda file_name: int(file_name.split('step')[0][10:]))
            last_state_of_mmp2_filepath = os.path.join(mmp2_files_path,mmp2_files[-1])
            last_state_of_ecm_filepath  = os.path.join(ecm_files_path,ecm_files[-1])
            print(f"Loading MMP2 state for grid id {grid_number + 1} in {last_state_of_mmp2_filepath}.")
            print(f"Loading ECM state for grid id {grid_number + 1} in {last_state_of_ecm_filepath}.")
            self.ecm[grid_number][0,:,:]  = pd.read_csv(last_state_of_ecm_filepath, index_col=0).to_numpy(dtype=float)
            self.mmp2[grid_number][0,:,:] = pd.read_csv(last_state_of_mmp2_filepath, index_col=0).to_numpy(dtype=float)

        path = os.path.join(path_to_simulation, "CellsData.csv")
        previous_sim_df = pd.read_csv(path, converters={"Position": ast.literal_eval})
        last_step = previous_sim_df["Step"].max()
        self.loaded_max_step = last_step
        previous_sim_df = previous_sim_df[previous_sim_df["Step"] == last_step]
        last_step_cells = previous_sim_df[previous_sim_df["Agent Type"] == "cell"]
        last_step_vessels = previous_sim_df[previous_sim_df["Agent Type"] == "vessel"]
        self.number_of_initial_cells = 0
        for index, row in last_step_cells.iterrows():
            current_grid_number = int(row["Grid"]) - 1
            ccell = CancerCell(self.current_agent_id, self, self.grids[current_grid_number], self.grid_ids[current_grid_number], row["Phenotype"], self.ecm[current_grid_number], self.mmp2[current_grid_number])
            self.current_agent_id += 1
            self.schedule.add(ccell)
            self.grids[current_grid_number].place_agent(ccell, row["Position"])
            self.cancer_cells_counter[current_grid_number] += 1
            self.number_of_initial_cells += 1
        for index, row in last_step_vessels.iterrows():
            current_grid_number = int(row["Grid"]) - 1
            ruptured_state = bool(row["Ruptured"])
            vessel = Vessel(self.current_agent_id, self, ruptured_state, self.grids[current_grid_number], self.grid_ids[current_grid_number])
            self.current_agent_id += 1
            self.schedule.add(vessel)
            self.grids[current_grid_number].place_agent(vessel, row["Position"])
            self.grid_vessels_positions[current_grid_number] += [row["Position"]]

        #load vasculature
        vasculature_path = os.path.join(path_to_simulation, "Vasculature")
        vasculature_files = os.listdir(vasculature_path)
        vasculature_files.sort(key = lambda file_name: int(file_name.split('step')[0][12:]))
        last_state_of_vasculature_filepath = os.path.join(vasculature_path,vasculature_files[-1])
        with open(last_state_of_vasculature_filepath, 'r') as f:
            last_state_of_vasculature = json.load(f)
        # Change keys to int
        last_state_of_vasculature = {int(k): v for k, v in last_state_of_vasculature.items()}
        self.vasculature = last_state_of_vasculature

        #calculate state of doubling counters
        self.doubling_time_counter_E = doubling_time_E - (last_step % doubling_time_E)
        self.doubling_time_counter_M = doubling_time_M - (last_step % doubling_time_M)

        #load time_grid_got_populated
        time_grid_got_populated_path = os.path.join(path_to_simulation, "Time when grids were populated")
        time_grid_got_populated_files = os.listdir(time_grid_got_populated_path)
        time_grid_got_populated_files.sort(key = lambda file_name: int(file_name.split('step')[0][25:]))
        time_grid_got_populated_filepath = os.path.join(time_grid_got_populated_path,time_grid_got_populated_files[-1])
        df_time_grid_got_populated = pd.read_csv(time_grid_got_populated_filepath, index_col=0)
        self.time_grid_got_populated = df_time_grid_got_populated.loc[0, :].values.flatten().tolist()


    def _initialize_grids(self):
        """
        Places the initial cancer cell and vessel in the initial grid in a circle

        Input: none
        Returns: none
        """
        mesenchymal_number = round(self.number_of_initial_cells * mesenchymal_proportion)
        possible_places = find_quasi_circle(n_center_points_for_tumor, self.width, self.height)[1]
        # Place all the agents in the quasi-circle area in the center of the grid
        for i in range(self.number_of_initial_cells):
            if mesenchymal_number > 0:
                cell_type = "mesenchymal"
                mesenchymal_number -= 1
            elif mesenchymal_number == 0:
                cell_type = "epithelial"

            a = CancerCell(self.current_agent_id, self, self.grids[0], self.grid_ids[0], cell_type, self.ecm[0], self.mmp2[0])
            self.current_agent_id += 1
            j = self.random.randrange(len(possible_places))
            x = int(possible_places[j][0])
            y = int(possible_places[j][1])

            self.schedule.add(a)
            self.grids[0].place_agent(a, (x, y))
            self.cancer_cells_counter[0] += 1

            # Remove the point after it has an amount of cells equal to the carrying capacity
            possible_places[j][2] += 1
            if possible_places[j][2] == carrying_capacity:
                possible_places.pop(j)


        # Create agents at second grid. Useful for debugging.
        amount_of_second_grid_cancer_cells = 0
        for i in range(amount_of_second_grid_cancer_cells):
            a = CancerCell(self.current_agent_id, self, self.grids[1], self.grid_ids[1], "mesenchymal", self.ecm[1], self.mmp2[1])
            self.current_agent_id += 1
            self.schedule.add(a)
        
            # Add the agent to a random grid cell
            x = self.random.randrange(3,7)
            y = self.random.randrange(3,7)
            self.grids[1].place_agent(a, (x, y))
            self.cancer_cells_counter[1] += 1

        # Create vessels
        num_normal_vessels = normal_vessels_primary
        num_ruptured_vessels = ruptured_vessels_primary

        # creates grid with 1 where vessels must not be placed
        not_possible_array = find_quasi_circle(n_center_points_for_Vessels, self.width, self.height)[0]
        not_possible_array[:2,:] = 1
        not_possible_array[-2:,:] = 1
        not_possible_array[:,:2] = 1
        not_possible_array[:,-2:] = 1
        possible_places = np.where(not_possible_array == 0)
        pos_coords = [list(tup) for tup in zip(possible_places[0], possible_places[1])]

        for i in range(len(self.grids)):
            if i == 0: # primary grid
                temp = num_ruptured_vessels
                while temp > 0:
                    j = num_ruptured_vessels - temp
                    coord_to_place = [self.random.randrange(self.width), self.random.randrange(self.height)]
                    if coord_to_place in pos_coords:
                        a = Vessel(self.current_agent_id, self, True, self.grids[0], self.grid_ids[0])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        self.grids[0].place_agent(a, (int(coord_to_place[0]), int(coord_to_place[1])))
                        self.grid_vessels_positions[i] += [(int(coord_to_place[0]), int(coord_to_place[1]))]
                        not_possible_array[coord_to_place[0], coord_to_place[1]] = 1
                        pos_coords.remove(coord_to_place)
                        temp -= 1

                temp = num_normal_vessels
                while temp > 0:
                    j = num_normal_vessels - temp
                    coord_to_place = [self.random.randrange(self.width), self.random.randrange(self.height)]
                    if coord_to_place in pos_coords:
                        a = Vessel(self.current_agent_id, self, False, self.grids[0], self.grid_ids[0])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        self.grids[0].place_agent(a, (int(coord_to_place[0]), int(coord_to_place[1])))
                        self.grid_vessels_positions[i] += [(int(coord_to_place[0]), int(coord_to_place[1]))]

                        not_possible_array[coord_to_place[0], coord_to_place[1]] = 1
                        pos_coords.remove(coord_to_place)
                        temp -= 1
            elif i > 0: # secondary grid and beyond
                    for m in range(secondary_sites_vessels[i-1]):
                        a = Vessel(self.current_agent_id, self, False, self.grids[i], self.grid_ids[i])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        x = self.random.randrange(self.width)
                        y = self.random.randrange(self.height)
                        self.grids[i].place_agent(a, (x,y))
                        self.grid_vessels_positions[i] += [(x,y)]
                
    def calculate_environment(self, mmp2, ecm):
        global th
        for i in range(len(mmp2)):
            for cell in self.grids[i].coord_iter():
                cell_contents, (x, y) = cell
                self.mesenchymal_count[i][x,y] = 0
                self.epithelial_count[i][x,y] = 0
                for cancer_cell in cell_contents:
                    if isinstance(cancer_cell, CancerCell):
                        if cancer_cell.phenotype == "mesenchymal":
                            self.mesenchymal_count[i][x,y] += 1
                        elif cancer_cell.phenotype == "epithelial":
                            self.epithelial_count[i][x,y] += 1
                        else:
                            raise Exception("Unknown phenotype")
                on_left_border = self.grids[i].out_of_bounds((x-1,y))
                on_right_border = self.grids[i].out_of_bounds((x+1,y))
                on_top_border = self.grids[i].out_of_bounds((x,y-1))
                on_bottom_border = self.grids[i].out_of_bounds((x,y+1))
                mmp2[i][1,x,y]=dmmp*tha/xha**2*((mmp2[i][0,x+1,y] if not on_right_border else mmp2[i][0,x-1,y])\
                        +(mmp2[i][0,x-1,y] if not on_left_border else mmp2[i][0,x+1,y])\
                        +(mmp2[i][0,x,y+1] if not on_bottom_border else mmp2[i][0,x,y-1])\
                        +(mmp2[i][0,x,y-1] if not on_top_border else mmp2[i][0,x,y+1])\
                        )\
                        +mmp2[i][0,x,y]*(1-4*dmmp*tha/xha**2-th*Lambda)+tha*theta*self.mesenchymal_count[i][x,y]
                ecm[i][1,x,y] = ecm[i][0,x,y]*(1-tha*(gamma1*self.mesenchymal_count[i][x,y]+gamma2*mmp2[i][1,x,y]))
                if ecm[i][1,x,y] < 0:
                    warnings.warn(f"<0 ecm in [i][1,{x},{y}] is {ecm[i][1,x,y]}")
                if ecm[i][1,x,y] > 1:
                    warnings.warn(f">1 ecm in [i][1,{x},{y}] is {ecm[i][1,x,y]}")
                    print("ECM is greater than 1! Your MMP2 diffusion rate is probably too high")
            mmp2[i][0,:,:] = mmp2[i][1,:,:]
            ecm[i][0,:,:] = ecm[i][1,:,:]

    def disaggregate_clusters(self, time):
        """
        Dissagregates cells from clusters into single-cell clusters, according to
        the dissagregation probability

        Input:
            time: given time for which the dissagreggation will occur
        Returns: 
            None
        """
        big_clusters = [cluster for cluster in self.vasculature[time] if sum(cluster) > 1]
        new_vasculature = [cluster for cluster in self.vasculature[time] if sum(cluster) == 1]
        for cluster in big_clusters:
            new_mesenchymal, new_epithelial = cluster
            for ccell_type, ccells_amount in enumerate(cluster):
                for i in range(ccells_amount):
                    if self.random.random() > dissagreggation_prob:
                        if ccell_type == 0:
                            new_vasculature += [(1, 0)]
                            new_mesenchymal -= 1
                        if ccell_type == 1:
                            new_vasculature += [(0, 1)]
                            new_epithelial -= 1
            if new_mesenchymal + new_epithelial > 0:
                new_vasculature += [(new_mesenchymal,new_epithelial)]
        self.vasculature[time] = new_vasculature