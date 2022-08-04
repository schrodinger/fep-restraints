#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
from tasks.io_trajectory import load_trajectory, write_frames

if __name__ == "__main__":
    """ 
    Read selected frames from a trajectory.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='+', dest='input_file',  type=str, help='Name of the input file(s) in csv format determining topology+trajectory files.')
    parser.add_argument('-n', dest='frame_number', type=int, help='frame number (zero-based)', default=0)
    parser.add_argument('-o', dest='output_dir', type=str, help='name of the output dir', default='start_frames')        
    args = parser.parse_args()

    # TODO: Set up the list with the frame number for each trajectory    

    # Read the input files
    print("\n* - Reading the input files. - *\n")
    simulations = []
    for csv in args.input_file:
        simulations.append( pd.read_csv(csv) )
    simulations = pd.concat(simulations)
    print(simulations) 

    # Create the output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Get a list of all systems
    systems = np.unique(simulations['System_Name'])
    print(systems)

    for sys in systems:
        sim_sys = simulations[simulations['System_Name']==sys]
        sim_num = 0
        for cms_file, trj_file in zip(sim_sys['Topology'], sim_sys['Trajectory']):
            sim_num += 1
            msys_model, cms_model = topo.read_cms(cms_file)
            trj = traj.read_traj(trj_file) 
            fname = "%s-sim%02i-frame%03i"%(sys, sim_num, args.frame_number)
            write_frames(
                cms_model, trj, [args.frame_number],
                args.output_dir, frame_names=[fname]
                )
