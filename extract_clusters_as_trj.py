#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import pandas as pd
import numpy as np
import os.path

from trajectory_io import load_trajectory, write_frames, get_an_aid
      

if __name__ == "__main__":
    """ 
    Read selected frames from a trajectory.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='cms_file', type=str, help='cms file')
    parser.add_argument('-t', nargs='+', dest='trj_files', type=str, help='trajecotry files or directories')    
    parser.add_argument('-n', nargs='+', dest='frame_number_files', type=str, help='files with frame numbers (zero-based) and properties')
    parser.add_argument('-p', dest='property', type=str, help='name of the column where the cluster ID is saved', default='Cluster_ID')
    parser.add_argument('-o', dest='output_name', type=str, help='name of the output files')  
    args = parser.parse_args()

    # Find all cluster IDs
    cl_ids = []
    for csv in args.frame_number_files:
        num_df = pd.read_csv(csv)
        cl_ids.append(num_df[args.property])
    cl_ids = np.unique(cl_ids)

    # Get the model and save it
    msys_model, cms_model = topo.read_cms(args.cms_file)
    #aid_list = cms_model.select_atom(args.selection)
    with structure.StructureWriter(args.output_name+'.cms') as writer:
        writer.append(cms_model.fsys_ct)

    # Go through all cluster IDs
    for value in cl_ids:
        # Find the numbers of the frames in the cluster
        frame_numbers = []
        for csv in args.frame_number_files:
            frame_num_csv = []
            num_df = pd.read_csv(csv)
            clust_id = num_df[args.property]
            frame_id = np.array(num_df.index)
            for fr, cl in zip(frame_id, clust_id):
                if cl == value:
                    frame_num_csv.append(fr)
            frame_numbers.append(frame_num_csv)
        # Find and write the frames corresponding to the selected numbers 
        frame_list = []
        for trj, num in zip( args.trj_files, frame_numbers ):
            trajectory = traj.read_traj(trj)         
            for f, frame in enumerate(trajectory):
                # Only proceed if the current frame is among the selected
                if f in num:
                    frame_list.append(frame)
        traj_name = args.output_name+'_cluster%02i'%value+'.xtc'
        traj.write_traj(frame_list, traj_name)
