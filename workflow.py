#!/usr/bin/env python3

import pandas as pd
from schrodinger.application.desmond.packages import analysis
from io_trajectory import load_trajectory, get_an_aid
from io_features import write_features_to_csv
from read_ca_distances_from_trajectory import calculate_ca_distances


if __name__ == "__main__":
    """ 
    Select cluster centroids, averages, and RMSF via PCA and k-means clustering on C-alpha distances between all selected residues. 
    """
       
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_file', type=str, help='Name of the input file in csv format determining topology+trajectory files.')
    parser.add_argument('-s', dest='selection_file', type=str, help='Name of the input file in csv format determining residue selections for binding pocket and restraints.') 
    parser.add_argument('-o', dest='output_dir', type=str, help='Name of the output directory', default='clustering_results')            
    args = parser.parse_args()

    trajectories = pd.read_csv(args.input_file)
    selections = pd.read_csv(args.selection_file)

    print(trajectories)
    print(selections)