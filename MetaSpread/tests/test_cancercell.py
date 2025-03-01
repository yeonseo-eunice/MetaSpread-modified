from metaspread.cancermodel import CancerModel
from metaspread.cancercell import CancerCell
import numpy as np
import pandas as pd
import pytest
import ast

#todo: model is not callable (duh! I think I cannot call a private variable (is it though?))
#todo: use tmp_path_facorty to create the model once, and use it for the rest of the tests

def test_cancercell(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=30,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder
        )
    
    assert model.data_collection_period == 10
    assert model.number_of_initial_cells==30
    assert model.width==201
    assert model.height==201
    assert model.grids_number==2
    assert model.max_steps==1000
    assert model.data_collection_period==10
    assert model.new_simulation_folder==temp_simulation_folder

    ccell_id = None
    grid_id = 1
    #add to the following line the keyword arguments for the cancer cell
    ccell = CancerCell(
        unique_id=ccell_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )

    assert ccell.unique_id          == ccell_id
    assert ccell.model              == model
    assert ccell.grid               == model.grids[grid_id-1]
    assert ccell.grid_id            == grid_id
    assert ccell.phenotype          == "mesenchymal"
    assert np.array_equal(ccell.ecm, model.ecm[grid_id-1])
    assert np.array_equal(ccell.mmp2, model.mmp2[grid_id-1])
def test_cancercell_movement_left(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=1,
        fixed_p_right=0,
        fixed_p_top=0,
        fixed_p_bottom=0
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        future_positions = [(pos[0]-1,pos[1]) for pos in current_positions]
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == future_positions[i]
                i = i + 1

def test_cancercell_movement_right(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=0,
        fixed_p_right=1,
        fixed_p_top=0,
        fixed_p_bottom=0
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        future_positions = [(pos[0]+1,pos[1]) for pos in current_positions]
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == future_positions[i]
                i = i + 1

def test_cancercell_movement_down(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=0,
        fixed_p_right=0,
        fixed_p_top=0,
        fixed_p_bottom=1
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        future_positions = [(pos[0],pos[1]-1) for pos in current_positions]
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == future_positions[i]
                i = i + 1

def test_cancercell_movement_up(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=0,
        fixed_p_right=0,
        fixed_p_top=1,
        fixed_p_bottom=0
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        future_positions = [(pos[0],pos[1]+1) for pos in current_positions]
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == future_positions[i]
                i = i + 1

def test_cancercell_movement_stay(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=0,
        fixed_p_right=0,
        fixed_p_top=0,
        fixed_p_bottom=0
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == current_positions[i]
                i = i + 1