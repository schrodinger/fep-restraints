#!/usr/bin/env python3

from utils.io_trajectory import load_trajectory
from utils.io_features import write_features_to_csv, calculate_ca_distances


if __name__ == "__main__":
    """ 
    Calculates the C-alpha distancses between all selected residues. 
    """
       
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='cms_file', type=str, help='cms file')
    parser.add_argument('-x', dest='xtc_file', type=str, help='xtc file', default=None)
    parser.add_argument('-t', dest='trj_dir', type=str, help='trajecotry directory', default=None)    
    parser.add_argument('-n', nargs='+', dest='residue_names', type=str, help='residue names', default=None)
    parser.add_argument('-r', nargs='+', dest='residue_numbers', type=int, help='residue numbers', default=None)
    parser.add_argument('-l', dest='chain_id', type=str, help='chain identifier', default='A')
    parser.add_argument('-o', dest='outputfile', type=str, help='name of the output csv file', default='results.csv')            
    args = parser.parse_args()

    # load data
    msys_model, cms_model, trajectory = load_trajectory(
        args.cms_file, 
        args.xtc_file, 
        args.trj_dir
        )
    # calculate the distances between the C-alpha atoms
    time, dist_names, distances = calculate_ca_distances(
        msys_model, cms_model, trajectory, 
        args.chain_id, args.residue_numbers, 
        residue_names=args.residue_names
        )
    # write the distances to a CSV file
    write_features_to_csv(
        args.outputfile, 
        distances, 
        dist_names, 
        time
        )
