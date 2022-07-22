#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
from io_trajectory import load_trajectory, write_frames

if __name__ == "__main__":
    """ 
    Read selected frames from a trajectory.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='cms_file', type=str, help='cms file')
    parser.add_argument('-x', dest='xtc_file', type=str, help='xtc file', default=None)
    parser.add_argument('-t', dest='trj_dir', type=str, help='trajecotry directory', default=None)    
    parser.add_argument('-n', nargs='+', dest='frame_numbers', type=int, help='frame numbers (zero-based)', default=[5,7])
    parser.add_argument('-o', dest='output_dir', type=str, help='name of the output dir', default='.')        
    args = parser.parse_args()

    msys_model, cms_model, trj = load_trajectory(args.cms_file, args.xtc_file, args.trj_dir)
    write_frames(cms_model, trj, args.frame_numbers, args.output_dir)
