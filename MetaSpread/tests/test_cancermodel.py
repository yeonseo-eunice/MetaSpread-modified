import pytest
from metaspread import cancermodel
from metaspread import configs
from metaspread import simrunner
from pathlib import Path

def test_imports():
    assert cancermodel.mesa is not None
    assert cancermodel.plt is not None
    assert cancermodel.np is not None
    assert cancermodel.pd is not None
    assert cancermodel.os is not None
    assert cancermodel.json is not None
    assert cancermodel.ast is not None
    assert cancermodel.CancerCell is not None
    assert cancermodel.Vessel is not None
    assert cancermodel.find_quasi_circle is not None

# def test_generate_cancer_model(mocker):
#     mocker.patch('metaspread.cancermodel.generate_cancer_model')
#     cancermodel.generate_cancer_model()
#     assert cancermodel.generate_cancer_model.called

from metaspread.cancermodel import CancerModel
from metaspread.cancercell import CancerCell
from metaspread.vessel import Vessel
import numpy as np
import pandas as pd
import pytest
import ast

#todo: model is not callable (duh! I think I cannot call a private variable (is it though?))
#todo: use tmp_path_facorty to create the model once, and use it for the rest of the tests

def test_load(tmp_path) -> None:
    simulations_path = tmp_path / "Simulations"
    model = simrunner.run_simulation(10, 10, save_path=tmp_path)
    
    model_loaded = CancerModel(
        number_of_initial_cells=30,
        width=201,
        height=201,
        grids_number=3,
        max_steps=1000,
        data_collection_period=2,
        new_simulation_folder= simulations_path / "test_loaded_sim",
        loaded_simulation_path = simulations_path / "Sim-max_steps-10-collection_period-10-cells-388-grids_number-3")
    assert model.vasculature == model_loaded.vasculature
    assert model.number_of_initial_cells == model_loaded.number_of_initial_cells
    assert model.width == model_loaded.width
    assert model.height == model_loaded.height
    assert model.phenotypes == model_loaded.phenotypes
    assert model.grid_vessels_positions == model_loaded.grid_vessels_positions
    assert model.current_agent_id == model_loaded.current_agent_id
    assert model.grids_number == model_loaded.grids_number
    assert model.grid_ids == model_loaded.grid_ids
    assert model.cancer_cells_counter == model_loaded.cancer_cells_counter
    assert model.time_grid_got_populated == model_loaded.time_grid_got_populated
    
    for i in range(1,model.grids_number):
        print(i)
        assert model.grid_vessels_positions[i] == model_loaded.grid_vessels_positions[i]
    
    for i in range(model.grids_number):
        assert (np.abs(model.mmp2[i][0] - model_loaded.mmp2[i][0]) < 10e-16).all()
        assert (np.abs(model.ecm[i][0] - model_loaded.ecm[i][0]) < 10e-16).all()
    
    assert model.max_steps != model_loaded.max_steps
    # right now the cells are counted on each step, so loading these cell counts
    # from the previous model is not implemented, therefore this test is not necessary.
    # But in the future that optimization should be implemented and thus,
    # this tests should be added 
    # for i in range(model.grids_number):
    #     assert (model.mesenchymal_count[i] == model_loaded.mesenchymal_count[i]).all()
    #     assert (model.epithelial_count[i] == model_loaded.epithelial_count[i]).all()


def test_cancermodel(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "Simulations" 
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=30,
        width=201,
        height=201,
        grids_number=3,
        max_steps=1000,
        data_collection_period=200000,
        new_simulation_folder=temp_simulation_folder)
    
    assert model.number_of_initial_cells==30
    assert model.width==configs.gridsize
    assert model.height==configs.gridsize
    assert model.grids_number==configs.grids_number
    assert model.max_steps==1000
    assert model.data_collection_period==200000
    assert model.new_simulation_folder==temp_simulation_folder

    #test proliferation
    cancermodel.carrying_capacity = 200
    # current_cell_count = list(map(type, model.schedule.agents)).count(CancerCell)
    current_cell_count = cancermodel.count_total_cells(model)
    model.proliferate("mesenchymal")
    model.proliferate("epithelial")
    assert 2*current_cell_count == list(map(type, model.schedule.agents)).count(CancerCell)
    current_cell_count = 2*current_cell_count

    #test calculate_environment
    grids_number = model.grids_number
    for i in range(grids_number):
        assert (model.ecm[i]  == 1.0).all() == True
        assert (model.mmp2[i] == 0.0).all() == True
    model.calculate_environment(model.mmp2, model.ecm)
    assert (model.ecm[0]  == 1.0).all() == False
    assert (model.mmp2[0] == 0.0).all() == False
    for i in range(1,grids_number):
        assert (model.ecm[i]  == 1.0).all() == True 
        assert (model.mmp2[i] == 0.0).all() == True
    
    #test count_vasculature_cells
    assert cancermodel.count_vasculature_cells(model) == 0
    model.vasculature = {1: [(100,100)]}
    assert cancermodel.count_vasculature_cells(model) == 200

    #test cell disaggregation
    for i in range(100):
        for j in range(100):
            model.vasculature = {1: [(i,j)]}
            model.disaggregate_clusters(1)
            assert sum(x+y for x,y in model.vasculature[1]) == i + j

    #test cell travel
    model.vasculature = {1: [(100,100)]}
    cancermodel.single_cell_survival = 1
    cancermodel.cluster_survival = 1
    model.step()
    model.step()
    assert list(map(type, model.schedule.agents)).count(CancerCell) == current_cell_count + 200
    current_cell_count = current_cell_count + 200
    assert current_cell_count == list(map(type, model.schedule.agents)).count(CancerCell)
    

    
