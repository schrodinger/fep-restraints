#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import pandas as pd
import numpy as np
import os.path

from io_trajectory import load_trajectory, write_frames, get_an_aid
      

def extract_cluster_as_xtc(trj_files, output_name, csv_files, value, property='Cluster_ID'):
    frame_list = []
    for csv, trj in zip(csv_files, trj_files):
        num_df = pd.read_csv(csv)
        clust_id = num_df[property]
        trajectory = traj.read_traj(trj)
        print( 'Length of CSV file:', len(clust_id), ' Length of trajectory:', len(trajectory))
        assert len(clust_id) == len(trajectory)
        for clid, frame in zip(clust_id, trajectory):
            if clid == value:
                frame_list.append(frame)
    traj_name = output_name+'_cluster%02i'%value+'.xtc'
    traj.write_traj(frame_list, traj_name)


if __name__ == "__main__":
    """ 
    Read selected frames from a trajectory.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='cms_file', type=str, help='cms file')
    parser.add_argument('-t', nargs='+', dest='trj_files', type=str, help='trajectory files or directories')    
    parser.add_argument('-n', nargs='+', dest='frame_number_files', type=str, help='files with frame numbers (zero-based) and properties')
    parser.add_argument('-p', dest='property', type=str, help='name of the column where the cluster ID is saved', default='Cluster_ID')
    parser.add_argument('-d', dest='definitions_file', type=str, help='csv file with definitions of cluster IDs, input files, and centroids. Input trajectories must be passed in the order that corresponds to the feature files', default=None)
    parser.add_argument('-o', dest='output_name', type=str, help='name of the output files')  
    parser.add_argument('-s', dest='selection', type=str, help='select which atoms to write', default='protein and ligand')
    args = parser.parse_args()


    # Write all CMS files
    _, cms_model = topo.read_cms(args.cms_file)
    #aid_list = cms_model.select_atom(args.selection)
    model = cms_model.copy()
    with structure.StructureWriter(args.output_name+'.cms') as writer:
        writer.append(model.fsys_ct)
        for st in model.comp_ct:
            writer.append(st)

    # Find all cluster IDs
    if args.definitions_file is None:
        cl_ids = []
        for csv in args.frame_number_files:
            num_df = pd.read_csv(csv)
            cl_ids.append(num_df[args.property])
        cl_ids = np.unique(cl_ids)
    else:
        def_df = pd.read_csv(args.definitions_file)
        cl_ids = np.array(def_df['Cluster_ID']) 
    print('Cluster IDs:', cl_ids)

    # Go through all cluster IDs and extract the corresponding frames
    for value in cl_ids:
        extract_cluster_as_xtc(args.trj_files, args.output_name, args.frame_number_files, value, property=args.property)

