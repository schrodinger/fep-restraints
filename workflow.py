#!/usr/bin/env python3

import os, os.path
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from schrodinger.application.desmond.packages import analysis, traj, topo
from io_trajectory import load_trajectory, get_an_aid
from io_features import write_features_to_csv, read_features_from_csv_files
from read_ca_distances_from_trajectory import calculate_ca_distances
from clustering_on_pca import kmeans_on_pca


if __name__ == "__main__":
    """ 
    Select cluster centroids, averages, and RMSF via PCA and k-means clustering on C-alpha distances between all selected residues. 
    """
       
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_file',  type=str, help='Name of the input file in csv format determining topology+trajectory files.')
    parser.add_argument('-s', dest='select_file', type=str, help='Name of the input file in csv format determining residue selections for binding pocket and restraints.') 
    parser.add_argument('-o', dest='output_dir',  type=str, help='Name of the output directory', default='clustering_results')       
    parser.add_argument('-n', dest='n_components',type=int, help='Number of top principal components to consider', default=3)     
    args = parser.parse_args()

    # Read the input files
    simulations = pd.read_csv(args.input_file)
    selections = pd.read_csv(args.select_file)
    print(simulations)
    print(selections)

    # Create the output directories
    os.makedirs(args.output_dir, exist_ok=True)
    for subdir in ['1-distances', '2-clustering', '3-rmsf']:
        newdir = os.path.join(args.output_dir, subdir)
        os.makedirs(newdir, exist_ok=True)  

    # * Calculate the C-alpha distances. *
    dist_files = []
    for i in simulations.index:
        # Get names and info for this simulation
        top_file = simulations['Topology'][i]
        trj_file = simulations['Trajectory'][i]
        sys_name = simulations['System_Name'][i]
        sys_prop = selections[selections['System_Name']==sys_name]
        chain_id = sys_prop['BindingPocket_Chain'].item()
        res_nums = [int(rn) for rn in sys_prop['BindingPocket_ResNum'].item().split(' ')] 
        # Load the simulation data
        msys_model, cms_model = topo.read_cms(top_file)
        trj = traj.read_traj(trj_file) 
        # Calculate the distances between the C-alpha atoms
        time, dist_names, distances = calculate_ca_distances(msys_model, cms_model, trj, chain_id, res_nums)
        # ... and write them to a CSV file
        out_file = os.path.join(args.output_dir,'1-distances/ca-distances_%04i.csv'%i)
        dist_files.append(out_file)
        write_features_to_csv(out_file, distances, dist_names, time)
    simulations['CA-Dist_File'] = dist_files

    # * Principal Component Analysis *
    names, data, origin, orig_id = read_features_from_csv_files(dist_files)
    pca = PCA(n_components=args.n_components)
    pca.fit(data.T)
    pc = pca.components_
    print(pca.explained_variance_ratio_) 

    # * K-Means Clustering *



