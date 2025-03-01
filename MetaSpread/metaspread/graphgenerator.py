import ast
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import sys
import metaspread.configs

def get_equally_spaced_array(passed_array, number_of_elems):
    passed_array = np.array(passed_array)
    indexes = np.round(np.linspace(0, len(passed_array)-1, number_of_elems)).astype(int)
    return list(zip(indexes,passed_array[indexes]))

def plot_cancer(fig_counter, grid_id, step, real_time_at_step, simulation_path, tumor_images_path):
    plt.style.use("seaborn-v0_8-darkgrid")
    figure_path = os.path.join(tumor_images_path, f'Cells-grid{grid_id}-step{step} - Tumor size at {real_time_at_step/(3600*24):.2f} days.png')
    if os.path.isfile(figure_path):
        return
    
    all_coords_data_path = os.path.join(simulation_path, "Data analysis", "Tumor dynamics")
    current_coords_filename = [file for file in os.listdir(all_coords_data_path) if file.startswith(f"Cells-grid{grid_id}-step{step} - Tumor size at")][0]
    current_coords_path = os.path.join(all_coords_data_path, current_coords_filename)
    coords_list = pd.read_csv(current_coords_path, index_col=0)
    Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr = coords_list.iloc[0], coords_list.iloc[1], coords_list.iloc[2], coords_list.iloc[3], coords_list.iloc[4], coords_list.iloc[5], coords_list.iloc[6], coords_list.iloc[7]
    plt.figure(fig_counter, figsize=(6, 6), facecolor='white')
    
    plt.scatter(Xm, Ym, marker='o', color='blue', alpha=0.5/metaspread.configs.carrying_capacity, label="Mesenchymal cells")
    plt.scatter(Xe, Ye, marker='h', color='orange', alpha=0.5/metaspread.configs.carrying_capacity, label="Epithelial cells")
    plt.scatter(Xv, Yv, marker='.', color='red', alpha=0.8, label="Vasculature points")
    plt.scatter(Xvr, Yvr, marker='+', color='darkred', alpha=0.8, label="Ruptured vasculature points")
    plt.xlim(0, metaspread.configs.gridsize)
    plt.ylim(0, metaspread.configs.gridsize)

    xticks = np.arange(0, metaspread.configs.gridsize, step=int(metaspread.configs.gridsize/6)) # 6 ticks
    xticklabels = [str(round(j,1)) for j in np.arange(0, 2.1, step = 2/201*(metaspread.configs.gridsize/6))]
    plt.xticks(xticks, xticklabels)
    plt.yticks(xticks, xticklabels)
    plt.xlabel("mm")
    plt.ylabel("mm")
    plt.grid(False)

    if step <= 2000:
        plt.legend(loc="upper left")
    plt.title(f'Tumor size at {real_time_at_step/(3600*24):.2f} days ({step} steps) - grid {grid_id}', fontsize = 13)

    # save the figure
    plt.savefig(figure_path)


def plot_growth_data(simulation_path, cells_images_path, grid_id, step, real_time_at_step):
    plt.style.use("seaborn-v0_8-darkgrid")
    path_to_save = os.path.join(cells_images_path, f'CellsGrowth-grid{grid_id}-step{step} - {real_time_at_step/(3600*24):.2f} days.png')
    if os.path.isfile(path_to_save):
        return
    all_cells_growth_data_path = os.path.join(simulation_path, "Data analysis", "Cells growth")
    current_growth_filename = [file for file in os.listdir(all_cells_growth_data_path) if file.startswith(f"CellsGrowth-grid{grid_id}-step{step} - ")][0]
    current_growth_path = os.path.join(all_cells_growth_data_path, current_growth_filename)
    df_cells_number = pd.read_csv(current_growth_path, index_col=0)

    plt.plot(df_cells_number["Days"], df_cells_number["Number of Mesenchymal Cells"], label="Mesenchymal cells", color='tab:blue')
    plt.plot(df_cells_number["Days"], df_cells_number["Number of Epithelial Cells"], label="Epithelial cells", color='tab:orange')

    plt.title(f'Cells growth - grid {grid_id}', fontsize = 13)
    plt.xlabel('Days')
    plt.ylabel('Number of cells')
    plt.legend(loc="upper left")
    plt.ylim(0)

    # save the figure
    plt.savefig(path_to_save)

def plot_MMP2_or_ECM(i, step, real_time_at_step, files_path, fig_counter, grid_id, path_to_save, type="Mmp2"):
    if type=="Mmp2":
        figure_path = os.path.join(path_to_save, f'{type}-grid{grid_id}-step{step} - {real_time_at_step/(3600*24):.2f} days.png')
    elif type=="Ecm":
        figure_path = os.path.join(path_to_save, f'{type}-grid{grid_id}-step{step} - {real_time_at_step/(3600*24):.2f} days.png')
    if os.path.isfile(figure_path):
        return
    try:
        df = pd.read_csv(files_path[i], index_col=0)
    except:
        raise Exception(f"Corrupted or empty file {files_path[i]}")
    plt.figure(fig_counter, figsize=(6, 5), facecolor='white')

    if type=="Mmp2":
        plt.imshow(df.T, vmin=0, vmax=3)
    elif type == "Ecm":
        plt.imshow(df.T , vmin=0, vmax=1)
        
    plt.colorbar()

    plt.xlim(0, metaspread.configs.gridsize)
    plt.ylim(0, metaspread.configs.gridsize)

    xticks = np.arange(0, metaspread.configs.gridsize, step=int(metaspread.configs.gridsize/6)) # 6 ticks
    xticklabels = [str(round(j,1)) for j in np.arange(0, 2.1, step = 2/201*(metaspread.configs.gridsize/6))]
    plt.xticks(xticks, xticklabels)
    plt.yticks(xticks, xticklabels)
    plt.xlabel("mm")
    plt.ylabel("mm")
    plt.grid(visible=None)

    plt.title(f'{type} at {real_time_at_step/(3600*24):.2f} days ({step} steps) - grid {grid_id}', fontsize = 13)
    
    plt.savefig(figure_path)

def plot_histogram(histogram_csv_file_path, all_histogram_images_path, step, real_time_at_step, grid_id):
    plt.style.use("seaborn-v0_8-darkgrid")
    path_to_save = os.path.join(all_histogram_images_path, f"Cells-grid{grid_id}-step{step} - Histogram at {real_time_at_step/(3600*24):.2f} days.png")
    if os.path.isfile(path_to_save):
        return
    histogram = pd.read_csv(histogram_csv_file_path, header=0, index_col=0)
    if histogram.empty:
        return
    plt.xlabel('Number of cells per grid-point')
    plt.ylabel('\nNr. grid-points with given nr. of cells')
    plt.title(f'Positions histogram at {real_time_at_step/(3600*24):.2f} days ({step} steps) - grid {grid_id}', fontsize = 13)
    plt.bar(histogram['Bins'], histogram['Frequency'])
    plt.xticks(range(metaspread.configs.carrying_capacity + 1))
    plt.xlim([-1, metaspread.configs.carrying_capacity + 1])
    plt.savefig(path_to_save)

def plot_vasculature_graphs(vasculature_df, path_to_save, max_step, real_delta_time):
    plt.style.use("seaborn-v0_8-darkgrid")
    figure_path = os.path.join(path_to_save, f'Vasculature-cells-step{max_step}.png')
    if not os.path.isfile(figure_path):
        # Prepare the data for the bar chart
        mesenchymal_data = vasculature_df["Mesenchymal cells"]
        epithelial_data = vasculature_df["Epithelial cells"]
        time_steps = vasculature_df["Time"]
        plt.style.use("seaborn-v0_8-darkgrid")

        #second plot, clusters
        plt.figure(facecolor='white')
        
        # Plot the data for each category
        plt.step(time_steps, mesenchymal_data, where='mid', color='tab:blue', label='Mesenchymal Cells', linewidth=1.5)
        plt.step(time_steps, epithelial_data, where='mid', color='tab:orange', label='Epithelial Cells', linewidth=1.5)

        # Set the chart labels, title, and legend
        plt.xlabel('Day')
        plt.ylabel('Cell Count')
        plt.title('Vasculature Cell Counts')
        plt.xticks([0, max_step/2, max_step], [0, f"{(max_step/2)*real_delta_time/(3600*24):.2f}", f"{max_step*real_delta_time/(3600*24):.2f}"])
        
        plt.xlim([-0.5, max_step + 0.5])
        plt.legend(loc="upper left")

        # save the figure
        plt.savefig(figure_path)
        plt.close()

    figure_path = os.path.join(path_to_save, f'Vasculature-clusters-step{max_step}.png')
    if not os.path.isfile(figure_path):
        total_cluster_data = vasculature_df["Total clusters"]
        multicellular_cluster_data = vasculature_df["Multicellular clusters"]
        time_steps = vasculature_df["Time"]
        #second plot, clusters
        plt.figure(facecolor='white')
        
        # Plot the data for each category
        plt.step(time_steps, total_cluster_data, where='mid', color='darkred', label='Total clusters', linewidth=1.5)
        plt.step(time_steps, multicellular_cluster_data, where='mid', color='green', label='Multicellular clusters', linewidth=1.5)

        # Set the chart labels, title, and legend
        plt.xlabel('Day')
        plt.ylabel('Cluster Count')
        plt.title('Vasculature Cluster Counts')
        plt.xticks([0, max_step/2, max_step], [0, f"{(max_step/2)*real_delta_time/(3600*24):.2f}", f"{max_step*real_delta_time/(3600*24):.2f}"])
        
        plt.xlim([-0.5, max_step + 0.5])
        plt.legend(loc="upper left")

        # save the figure
        plt.savefig(figure_path)
        plt.close()

def plot_radius_diameter_history(df, path_to_save, max_step, real_delta_time):
    plt.style.use("seaborn-v0_8-darkgrid")
    figure_path = os.path.join(path_to_save, 'Radius and diameter for Grid 1.png')
    if not os.path.isfile(figure_path):
        # Prepare the data for the bar chart
        radius = df["Radius"]*metaspread.configs.xh*(0.001/0.005)*10  #getting distance in mm
        diameter = df["Diameter"]*metaspread.configs.xh*(0.001/0.005)*10 #getting distance in mm
        time_steps = df["Step"]
        plt.style.use("seaborn-v0_8-darkgrid")

        #second plot, clusters
        plt.figure(facecolor='white')
        
        # Plot the data for each category
        plt.step(time_steps, radius, where='mid', color='tab:blue', label='Radius', linewidth=1.5)
        plt.step(time_steps, diameter, where='mid', color='tab:orange', label='Diameter', linewidth=1.5)

        # Set the chart labels, title, and legend
        plt.xlabel('Day')
        plt.ylabel('Distance (mm)')
        plt.title('Radius and diameter for Grid 1')
        plt.xticks([0, max_step/2, max_step], [0, f"{(max_step/2)*real_delta_time/(3600*24):.2f}", f"{max_step*real_delta_time/(3600*24):.2f}"])
        
        plt.xlim([-0.5, max_step + 0.5])
        plt.legend(loc="upper left")

        # save the figure
        plt.savefig(figure_path)
        plt.close()

    
def generate_graphs(name_of_the_simulation, amount_of_pictures=0):
    simulations_dir = "Simulations"
    simulation_path = os.path.join(simulations_dir, name_of_the_simulation)
    configs_path = os.path.join(simulation_path, "configs.csv")
    metaspread.configs.load_simulation_configs_for_data_generation(configs_path)
    print(f'Analyzing data in the folder {simulation_path}\n')

    # Get the Ecm and Mmp2 data filenames 
    ecm_path = os.path.join(simulation_path, "Ecm")
    mmp2_path = os.path.join(simulation_path, "Mmp2")
    ecm_files_name = [f for f in sorted(os.listdir(ecm_path), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(ecm_path, f)) and f.endswith(".csv")]
    mmp2_files_name = [f for f in sorted(os.listdir(mmp2_path), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(mmp2_path, f)) and f.endswith(".csv")]

    # Get the vasculature data filename
    vasculature_path = os.path.join(simulation_path, "Vasculature")
    vasculature_files_name = [f for f in os.listdir(vasculature_path) if os.path.isfile(os.path.join(vasculature_path, f)) and f.endswith(".json")]

    all_cells_filename = os.path.join(simulation_path, "CellsData.csv")
    if os.path.isfile(all_cells_filename):
        print("Using cells data at:", all_cells_filename)
    else:
        print("No CellsData.csv cell data found in directory:", simulation_path)
        return

    if ecm_files_name:
        ecm_files_path = [os.path.join(ecm_path, p) for p in ecm_files_name]
        print("Using Ecm data in the folder:", ecm_path)
    else:
        print("No .csv Ecm data found in directory:", ecm_path)
        return

    if mmp2_files_name:
        mmp2_files_path = [os.path.join(mmp2_path, p) for p in mmp2_files_name]
        print("Using Mmp2 data in the folder:", mmp2_path)
    else:
        print("No .csv Mmp2 data found in directory:", mmp2_path)
        return

    if vasculature_files_name:
        print("Using vasculature data at:", vasculature_path)
    else:
        print("No .json vasculature data found in directory:", simulation_path)
        return


    #Path of the data
    data_folder = "Data analysis"
    data_path = os.path.join(simulation_path, data_folder)

    tumor_data_path = os.path.join(data_path, "Tumor dynamics")
    
    step_size = metaspread.configs.data_collection_period
    real_delta_time = 40 * metaspread.configs.th/0.001 #in seconds (the original ratio is 40 seconds/0.001 non-dimensional time)
    grids_number = metaspread.configs.grids_number
    configs_max_step = metaspread.configs.max_steps
    all_cells_dataframe = pd.read_csv(all_cells_filename, converters={"Position": ast.literal_eval})
    all_cells_dataframe = all_cells_dataframe[["Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]
    max_step = max(all_cells_dataframe["Step"])
    if configs_max_step >= max_step:
        print(f"Warning: the run for this simulation terminated early")
        print(f"Max step reached is {max_step} while {configs_max_step} was expected.")
    maximum_frames = int(max_step / step_size)
    if amount_of_pictures > maximum_frames:
        print(f"There are only {maximum_frames} frames available!")
        print("Setting up the amount of picture to this value.")
        amount_of_pictures = 0 #to print all the pictures
    
    # Path to save all the images:
    images_folder = "Graphical analysis"
    images_path = os.path.join(simulation_path, images_folder)

    tumor_images_path = os.path.join(images_path, "Tumor dynamics")
    cells_images_path = os.path.join(images_path, "Cells growth")
    ecm_images_path = os.path.join(images_path, "Ecm dynamics")
    mmp2_images_path = os.path.join(images_path, "Mmp2 dynamics")
    vasculature_images_path = os.path.join(images_path, "Vasculature dynamics")
    all_histogram_images_path = os.path.join(images_path, "Positions histogram")
    radius_diameter_images_path = os.path.join(images_path, "Radius and diameter")

    # Create folder for all the Graphical analysis

    os.makedirs(images_path, exist_ok = True)
    os.makedirs(tumor_images_path, exist_ok = True)
    os.makedirs(mmp2_images_path, exist_ok = True)
    os.makedirs(ecm_images_path, exist_ok = True)
    os.makedirs(cells_images_path, exist_ok = True)
    os.makedirs(vasculature_images_path, exist_ok = True)
    os.makedirs(all_histogram_images_path, exist_ok = True)
    os.makedirs(radius_diameter_images_path, exist_ok = True)

    print(f"\nSaving tumor images in the folder:", tumor_images_path)
    print("Saving Mmp2 images in the folder:", mmp2_images_path)
    print("Saving Ecm images in the folder:", ecm_images_path)
    print("Saving cells numbers images in the folder:", cells_images_path)
    print("Saving vasculature images in the folder:", vasculature_images_path)
    print("Saving histogram images in the folder:", all_histogram_images_path)
    print("Saving histogram image in the folder:", radius_diameter_images_path)

    
    if amount_of_pictures != 0:
        range_of_pictures = get_equally_spaced_array(range(step_size,max_step+1,step_size), amount_of_pictures)
    else:
        range_of_pictures = enumerate(range(step_size,max_step+1,step_size))
    fig_counter = 1
    for grid_id in range(1, grids_number+1):
        plt.style.use("default")
        plt.style.use("Solarize_Light2")
        print(f'\nGrid: {grid_id}')
    
        # Plot the Mmp2 graphs
        print(f'\tPlotting Mmp2 graphs...')
        for id, step in range_of_pictures:
            real_time_at_step = real_delta_time * step
            mmp2_files_path_this_grid = [path for path in mmp2_files_path if f"Mmp2-{grid_id}grid-" in path]
            plot_MMP2_or_ECM(id, step, real_time_at_step, mmp2_files_path_this_grid, fig_counter, grid_id, mmp2_images_path, type="Mmp2")
            plt.close()
            fig_counter += 1

        # Plot the Ecm graphs
        print(f'\tPlotting Ecm graphs...')
        for id, step in range_of_pictures:
            real_time_at_step = real_delta_time * step
            ecm_files_path_this_grid = [path for path in ecm_files_path if f"Ecm-{grid_id}grid-" in path]
            plot_MMP2_or_ECM(id, step, real_time_at_step, ecm_files_path_this_grid, fig_counter, grid_id, ecm_images_path, type="Ecm")
            plt.close()
            fig_counter += 1

        # Plot the cells graphs
        print(f'\tPlotting tumor graphs...')
        for id, step in range_of_pictures:
            real_time_at_step = real_delta_time * step
            plot_cancer(fig_counter, grid_id, step, real_time_at_step, simulation_path, tumor_images_path)
            plt.close()
            fig_counter += 1

        # Plot the histogram graphs
        print(f'\tPlotting histogram graphs...')
        csv_histogram_files_names = [f for f in os.listdir(tumor_data_path) if os.path.isfile(os.path.join(tumor_data_path, f)) and "Histogram" in f]
        for id, step in range_of_pictures:
            real_time_at_step = real_delta_time * step
            histogram_csv_file_name = csv_histogram_files_names[id]
            histogram_csv_file_path = os.path.join(tumor_data_path, histogram_csv_file_name)
            plot_histogram(histogram_csv_file_path, all_histogram_images_path, step, real_time_at_step, grid_id)
            plt.close()
            fig_counter += 1

        # Plot the growth of ephitelial and mesenchymal cells
        print(f'\tPlotting cells numbers graph...')
        for id, step in range_of_pictures:
            real_time_at_step = real_delta_time * step
            plot_growth_data(simulation_path, cells_images_path, grid_id, step, real_time_at_step)
            plt.close()
            fig_counter += 1

    # Plot the vasculature data
    print(f'Plotting vasculature...')
    for id, step in range_of_pictures:
        folder_path = os.path.join(simulation_path, "Data analysis", "Vasculature dynamics")
        try:
            file_path = os.path.join(folder_path, f"Vasculature-step{step}.csv")
            vasculature_data = pd.read_csv(file_path, header=0)
        except:
            print(f"Error while reading the vasculature data for step {step} in grid {grid_id}", file=sys.stderr)
            print("Did you run the 'Data analysis' in the postprocessing menu first?", file=sys.stderr)
            os._exit(1)
        plot_vasculature_graphs(vasculature_data, vasculature_images_path, step, real_delta_time)
        plt.close()
        fig_counter += 1
        
    #plotting the radius and diameter history graph
    print(f'Plotting radius and diameter history graph...')
    folder_path = os.path.join(simulation_path, "Data analysis", "Tumor dynamics")
    file_name = [file for file in os.listdir(folder_path) if file.startswith("Tumor radius and diameter history")][0]
    file_path = os.path.join(folder_path, file_name)
    radius_history_df = pd.read_csv(file_path)#, header=0)
    plot_radius_diameter_history(radius_history_df, radius_diameter_images_path, step, real_delta_time)
    plt.close()
    fig_counter += 1

    plt.show()
    plt.close()

    plt.show()
    plt.close()






