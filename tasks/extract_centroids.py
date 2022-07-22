#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import pandas as pd
import numpy as np
import os.path

from utils.io_trajectory import load_trajectory, write_frames, get_an_aid
      

if __name__ == "__main__":
    """ 
    Read selected frames from a trajectory.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', nargs='+', dest='cms_files', type=str, help='cms files')
    parser.add_argument('-t', nargs='+', dest='trj_files', type=str, help='trajecotry files or directories')    
    parser.add_argument('-d', dest='definitions_file', type=str, help='csv file with definitions of cluster IDs, input files, and centroids. Input trajectories must be passed in the order that corresponds to the feature files', default=None)
    parser.add_argument('-o', dest='output_name', type=str, help='name of the output files')  
    args = parser.parse_args()

    # Read definitions file 
    def_df = pd.read_csv(args.definitions_file)
    cl_ids = np.array(def_df['Cluster_ID']) 
    c_files = np.array(def_df['COrig_File_ID'])
    c_frames = np.array(def_df['Centroid_Original_Index'])

    # Get the trajectory and save the frame with the centroid
    for cl_id, c_file, c_frame in zip( cl_ids, c_files, c_frames ):
        _, top = topo.read_cms(args.cms_files[c_file])
        trj = traj.read_traj(args.trj_files[c_file])
        out_dir = os.path.dirname(args.output_name)
        out_fn = os.path.basename(args.output_name) + '_centroid%02i'%cl_id
        write_frames(top, trj, [c_frame], out_dir, frame_names=[out_fn])
