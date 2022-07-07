#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import pandas as pd
import numpy as np
import os.path

from io_trajectory import load_trajectory, write_frames, get_an_aid


def select_subset_model(cms_model, aid):
    aid_all = range(1,cms_model.atom_total+1)
    aid_del = set(aid_all) - set(aid)
    model = cms_model.copy()
    model.deleteAtoms(aid_del)
    print('The new model has %s atoms.'%model.atom_total)
    print('The new model has %s atoms.'%len(model.atom)) 
    return model


def write_frame(out_fname, old_model, frame, sigma=None):
    model = old_model.copy()
    if frame is not None:      
        for i, ct in enumerate(model.comp_ct):
            ct_aids = model.get_fullsys_ct_atom_index_range(i)
            ct_gids = model.convert_to_gids(ct_aids, with_pseudo=False)
            topo.update_ct(ct, model, frame, is_fullsystem=False, allaid_gids=ct_gids)
        model.synchronize_fsys_ct()
    # Add the RMSF as the atom sigma
    if sigma is not None:
        for a, atom in enumerate(model.atom):
            atom.property['r_desmond_sigma'] = sigma[a]
    # Write the structure to the designated name
    with structure.StructureWriter(out_fname) as writer:
        writer.append(model.fsys_ct)
        for st in model.comp_ct:
            writer.append(st)
    return model


def write_coordinates(out_fname, cms_model, xyz=None, sigma=None):
    model = cms_model.copy()  
    # Update the coordinates. 
    if xyz is not None:
        model.fsys_ct.setXYZ(xyz)
    # Add the RMSF as the atom sigma
    if sigma is not None:
        for a, atom in enumerate(model.atom):
            atom.property['r_desmond_sigma'] = sigma[a]
    # Write the structure to the designated name.
    with structure.StructureWriter(out_fname) as writer:
        writer.append(model.fsys_ct)


if __name__ == "__main__":
    """ 
    Read selected frames from a trajectory.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', nargs='+', dest='cms_files', type=str, help='cms files')
    parser.add_argument('-t', nargs='+', dest='trj_files', type=str, help='trajecotry files or directories')    
    parser.add_argument('-s', dest='sel_file', type=str, help='alignment selections for each trajectory, each line is one selection', default='selections.txt')
    parser.add_argument('-w', dest='write_sel_file', type=str, help='selection file for which atoms to use in each trajectory')
    parser.add_argument('-o', dest='output_filename', type=str, help='name of the output files')  
    parser.add_argument('--ref_file', dest='reference_fn', type=str, help='reference file to align and calculate RMSF to')
    parser.add_argument('--ref_sel', dest='reference_asl', type=str, help='alignment selection for the reference file')
    args = parser.parse_args()

    # Read the reference structure
    _, cms_model_ref = topo.read_cms(args.reference_fn) 
    aidlist_ref = cms_model_ref.select_atom(args.reference_asl) 
    gidlist_ref = topo.aids2gids(cms_model_ref, aidlist_ref, include_pseudoatoms=False)
    pos_ref = cms_model_ref.getXYZ()[gidlist_ref]

    # Load selection files
    with open(args.sel_file) as sf:
        selections_align = [line.rstrip() for line in sf]
    with open(args.write_sel_file) as sf:
        selections_write = [line.rstrip() for line in sf]

    # Create list in which to store the new positions and distances
    pos_sel_aligned = []
    squared_distances = []
    # Loop through all trajectories, align, and calculate distances
    for cms, trj, sel_align, sel_write in zip(args.cms_files, args.trj_files, selections_align, selections_write):
        # Read the trajectory
        _, cms_model = topo.read_cms(cms)
        print('Number of atoms in the system:', cms_model.atom_total)
        trajectory = traj.read_traj(trj)
        print('Alignment selection:', str(sel_align))
        aidlist_align = cms_model.select_atom(str(sel_align))
        gidlist_align = topo.aids2gids(cms_model, aidlist_align, include_pseudoatoms=False)
        print('Output Selection:', str(sel_write))
        aidlist_write = cms_model.select_atom(str(sel_write)) 
        gidlist_write = topo.aids2gids(cms_model, aidlist_write, include_pseudoatoms=False)
        # Loop through trajectory
        model = cms_model.copy()
        ct = model.fsys_ct
        print('Reading:', cms, trj)
        for f, frame in enumerate(trajectory):
            topo.update_ct(ct, model, frame)
            # Align frame to reference structure
            pos_all = ct.getXYZ()
            pos_sel = pos_all[gidlist_align]
            pos_new = analysis.align_pos(pos_all, pos_sel, pos_ref)
            # Calculate distances to reference structure and append them to the list
            pos_sel_aligned.append(pos_new[gidlist_write])
            squared_distances.append(np.mean((pos_new[gidlist_write]-pos_all[gidlist_write])**2,axis=1))
    pos_sel_aligned = np.array(pos_sel_aligned)
    squared_distances = np.array(squared_distances)

    # Calculate the average position of each selected atom
    pos_average = np.mean(pos_sel_aligned, axis=0)

    # Calculate the RMSF for each atom
    rmsf_per_atom = np.sqrt(np.mean(squared_distances, axis=0))
    rmsd_per_frame = np.sqrt(np.mean(squared_distances, axis=1))
    
    output = pd.DataFrame(rmsf_per_atom)
    output.to_csv(args.output_filename+'.csv')
    
    # Write RMSF on reference structure
    #out_fn_ref = args.output_filename+'_rmsf_ref.cms'
    #_ = write_coordinates(out_fn_ref, cms_model_ref, xyz=None, sigma=rmsf_per_atom)

    # Write RMSF on average structure
    #out_fn_avg = args.output_filename+'_rmsf_avg.cms'
    #_ = write_coordinates(out_fn_avg, cms_model_ref, pos_average, sigma=rmsf_per_atom)
    