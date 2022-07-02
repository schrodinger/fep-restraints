#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import pandas as pd
import os.path


def get_an_aid(cms_model, asl_string):
    aid_list = cms_model.select_atom(asl_string)
    if (len(aid_list) != 1):
        raise Exception("The number of chosen atoms is %d."%(len(aid_list)))
    else:
        return(aid_list[0])

def load_trajectory(cms_file, xtc_file, trj_dir):
    msys_model, cms_model = topo.read_cms(cms_file)
    if ((not xtc_file) and (not trj_dir)):
        raise Exception("Neither the xtc file nor the trajectory directory is found.")
    elif (xtc_file and trj_dir):
        raise Exception("Both the xtc file and the trajectory directory are found.")
    elif (xtc_file):
        trj = traj.read_traj(xtc_file)  # xtc format trajectory
    else:
        trj = traj.read_traj(trj_dir)  # trj format trajectory
    return cms_model, trj

def write_frames(cms_model, trajectory, frame_numbers, out_dir, frame_names=None):
    
    # Determine the names of the frames (without extension)
    if frame_names is None:
        frame_names = { n:'frame_%04i'%n for n in frame_numbers }
    else:
        assert len(frame_names) == len(frame_numbers)
        frame_names = { n: frame_names[i] for i,n in enumerate(frame_numbers)}
    # Work with a copy of the ct
    #fsys_ct = cms_model.fsys_ct.copy()
    model = cms_model.copy()
    # Iterate through the trajectory
    for f, frame in enumerate(trajectory):
        # Only proceed if the current frame is among the selected
        if f in frame_numbers:
            # Update the ct positions to the traj positions
            # We have to do this because traj includes pseudo-atoms
            for i, ct in enumerate(model.comp_ct):
                at_indices = model.get_fullsys_ct_atom_index_range(i)
                #at_gids = [ati-1 for ati in at_indices]
                at_gids = model.convert_to_gids(at_indices, with_pseudo=False)
                topo.update_ct(ct, model, frame, is_fullsystem=False, allaid_gids=at_gids)
            model.synchronize_fsys_ct()
            # Write the structure to the designated name
            out_fname = os.path.join(out_dir, frame_names[f])
            with structure.StructureWriter(out_fname+'.cms') as writer:
                writer.append(model.fsys_ct)
                for st in model.comp_ct:
                    writer.append(st)
    return



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

    cms_model, trj = load_trajectory(args.cms_file, args.xtc_file, args.trj_dir)
    write_frames(cms_model, trj, args.frame_numbers, args.output_dir)
