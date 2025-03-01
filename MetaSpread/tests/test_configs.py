import pytest
from metaspread import configs
from metaspread.cancermodel import CancerModel
import pandas as pd
import ast

def test_imports():
    assert configs.pd is not None
    assert configs.os is not None
    assert configs.ast is not None

def test_init_simulation_configs(mocker):
    mocker.patch('metaspread.configs.init_simulation_configs')
    configs.init_simulation_configs('mock_configs.csv')
    assert configs.init_simulation_configs.called


def test_configs(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation"
    temp_simulation_folder.mkdir()
    
    # Parameters for this simulation
    number_of_initial_cells = configs.number_of_initial_cells # Number of cancer cells
    gridsize     = configs.gridsize
    grids_number = configs.grids_number
    width        = gridsize
    height       = gridsize
    max_steps    = 1000
    data_collection_period = 10

    model = CancerModel(
        number_of_initial_cells=number_of_initial_cells,
        width=width,
        height=height,
        grids_number=grids_number,
        max_steps=max_steps,
        data_collection_period=data_collection_period,
        new_simulation_folder=temp_simulation_folder)
    
    #test against the csv file
    sim_configs = pd.read_csv('simulations_configs.csv', header=0, converters={"Values": ast.literal_eval})
    assert model.number_of_initial_cells== sim_configs.query("Names=='number_of_initial_cells'")["Values"].values[0]
    assert model.width == sim_configs.query("Names=='gridsize'")["Values"].values[0]
    assert model.height == sim_configs.query("Names=='gridsize'")["Values"].values[0]
    assert model.grids_number == sim_configs.query("Names=='grids_number'")["Values"].values[0]
    assert model.max_steps==max_steps
    assert model.data_collection_period == data_collection_period
    assert model.new_simulation_folder==temp_simulation_folder
