#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import pandas as pd
import numpy as np
import os.path

from io_trajectory import load_trajectory, write_frames, get_an_aid
      

if __name__ == "__main__":
    """ 
    Read selected frames from a trajectory.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', nargs='+', dest='cms_files', type=str, help='cms files')
    parser.add_argument('-t', nargs='+', dest='trj_files', type=str, help='trajecotry files or directories')    
    parser.add_argument('-s', nargs='+', dest='selection', type=str, help='asl strings for each trajectory')
    parser.add_argument('-o', dest='output_filename', type=str, help='name of the output files')  
    parser.add_argument('--ref_file', dest='reference_fn', type=str, help='reference file to align and calculate RMSF to')
    parser.add_argument('--ref_sel', dest='reference_asl', type=str, help='selection for the reference file')
    args = parser.parse_args()


    _, cms_model_ref = topo.read_cms(args.reference_fn) 
    aid_list_ref = cms_model_ref.select_atom(args.reference_asl) 
    coords_ref = cms_model_ref.getXYZ()[aid_list_ref]
    print('Reference coordinate:', coords_ref)

    # Create list in which to store the distances
    distances = []

    # Loop through all trajectories
    for cms, trj, sel in zip(args.cms_files, args.trj_files, args.selection):

        # Read the trajectory
        _, cms_model = topo.read_cms(cms)
        trajectory = traj.read_traj(trj)
        aid_list = cms_model.select_atom(sel)

        # Loop through trajectory
        model = cms_model.copy()
        for f, frame in enumerate(trajectory):
            print('Frame:', frame)
            topo.update_ct(model.fsys_ct, model, frame)

            # TODO: Align frame to reference structure
            # TODO: Calculate distances to reference structure and append them to the list

    # TODO: suare the distances for each atom
    # TODO: Calculate the average position of each atom
    # TODO: Calculate the RMSF for each atom

    # TODO: Write RMSF on average structure
    # TODO: Write RMSF on reference structure
     
    with structure.StructureWriter(args.output_name+'.cms') as writer:
        writer.append(cms_model.fsys_ct)

