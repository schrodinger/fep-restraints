#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .io_trajectory import load_trajectory, write_frames, get_an_aid, select_subset_model



def write_frame(out_fname, old_model, frame, sigma=None):
    model = old_model.copy()
    print("Writing to model with %i atoms."%(len(model.atom)))
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


def write_coordinates(out_fname, model, xyz=None, sigma=None):
    print("Writing to model with %i atoms."%(len(model.atom))) 
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


def calculate_rmsf(reference_fn, cms_files, trj_files, csv_files, cluster_id, ref_asl_align, ref_asl_write, selections_align, selections_write, 
                   align_avg=True, threshold=None, start_frame=0, end_frame=None, step=1):

    # Read the reference structure
    _, cms_model_ref = topo.read_cms(reference_fn) 
    aidlist_ref = cms_model_ref.select_atom(ref_asl_align) 
    gidlist_ref = topo.aids2gids(cms_model_ref, aidlist_ref, include_pseudoatoms=False)
    print("The reference selection has %i atoms."%(len(gidlist_ref)))
    pos_ref_all = cms_model_ref.getXYZ()
    pos_ref = pos_ref_all[gidlist_ref]

    # Create list in which to store the new positions and distances
    pos_sel_aligned = []
    squared_distances = []
    # Loop through all trajectories, align, and calculate distances
    for cms, trj, csv, sel_align, sel_write in zip(cms_files, trj_files, csv_files, selections_align, selections_write):
        # Read the CSV file with cluster information
        frame_cluster_id = np.array(pd.read_csv(csv)['Cluster_ID'], dtype=int)
        # Read the trajectory
        _, cms_model = topo.read_cms(cms)
        print('Number of atoms in the system:', cms_model.atom_total)
        trajectory = traj.read_traj(trj)[start_frame:end_frame:step]
        print('Alignment selection:', str(sel_align))
        aidlist_align = cms_model.select_atom(str(sel_align))
        gidlist_align = topo.aids2gids(cms_model, aidlist_align, include_pseudoatoms=False)
        print("Selected %i atoms for alignment."%(len(gidlist_align)))
        print('Output Selection:', str(sel_write))
        aidlist_write = cms_model.select_atom(str(sel_write)) 
        gidlist_write = topo.aids2gids(cms_model, aidlist_write, include_pseudoatoms=False)
        print("Selected %i atoms for output."%(len(gidlist_write)))
        # Loop through trajectory
        model = cms_model.copy()
        ct = model.fsys_ct
        print('Reading:', cms, trj)
        for f, frame in enumerate(trajectory):
            # Check whether this frame is in the selected cluster
            if frame_cluster_id[f] != cluster_id:
                continue
            # Update the coordinates
            topo.update_ct(ct, model, frame)
            # Align frame to reference structure
            pos_all = ct.getXYZ()
            pos_sel = pos_all[gidlist_align]
            pos_new = analysis.align_pos(pos_all, pos_sel, pos_ref)
            # Calculate distances to reference structure and append them to the list
            new_squ_dist = np.sum((pos_new[gidlist_write]-pos_ref_all[gidlist_write])**2, axis=1)
            frame_rmsd = np.sqrt(np.mean(new_squ_dist))
            if f%100 == 0: 
                print("RMSD in frame %04i: %1.4f"%(f, frame_rmsd))
            if threshold is None or frame_rmsd < threshold: 
                squared_distances.append(new_squ_dist)
                pos_sel_aligned.append(pos_new[gidlist_write])
    pos_sel_aligned = np.array(pos_sel_aligned)
    squared_distances_from_reference = np.array(squared_distances)
    print('Squared distances from reference:', squared_distances_from_reference.shape)

    # Calculate the average position of each selected atom.
    pos_average = np.mean(pos_sel_aligned, axis=0)
    # Calculate the RMSD of each frame to the average.
    squared_distances_from_average = np.sum((pos_sel_aligned-pos_average)**2, axis=2)
    print('Squared distances from average:', squared_distances_from_average.shape)
    # ... and align it to the reference, using all atoms
    if align_avg:
        pos_average = analysis.align_pos(pos_average, pos_average, pos_ref_all[gidlist_write])

    # Calculate the RMSF for each atom.
    rmsf_per_atom = np.sqrt(np.mean(squared_distances_from_average, axis=0))
    rmsd_per_frame = np.sqrt(np.mean(squared_distances_from_reference, axis=1))
    print("Calculated RMSF for %i atoms and %i frames."%(len(rmsf_per_atom), len(rmsd_per_frame)))    

    # Get the model of the subset to write to the structures.
    aidlist_write_ref = cms_model_ref.select_atom(str(ref_asl_write))
    cms_model_ref_new = select_subset_model(cms_model_ref, aidlist_write_ref )
    print('The new model has %s atoms.'%cms_model_ref_new.atom_total)
    print('The new model has %s atoms.'%len(cms_model_ref_new.atom)) 

    return rmsf_per_atom, pos_average, cms_model_ref_new


def plot_cluster_rmsf(k, rmsf_files_k, features, out_plot, highlight_feature_residues=True):
    # Plot the RMSF of each cluster
    fig, ax = plt.subplots(k, 1, figsize=[8,2.0*k], dpi=300, sharex=True, sharey=True)
    if k==1: # Make ax a list, even when there is only one plot
        ax = [ax]
    resnum_min = []
    resnum_max = []
    for cluster_id in range(k):
        # Get the data
        csv = rmsf_files_k[cluster_id]
        df = pd.read_csv(csv)
        df_ca = df[df['pdbname']==' CA ']
        # Find the range of residue numbers in this simulation
        resnum_min.append( np.min(df_ca['resnum']) )
        resnum_max.append( np.max(df_ca['resnum']) )
        # Draw a grid
        ax[cluster_id].grid(axis='y', ls='--')
        # Plot the RMSF of all residues
        ax[cluster_id].bar(
            df_ca['resnum'], df_ca['RMSF'], 
            label=r'C$\mathrm{\alpha}$ atoms', 
            width=1.0, color="C%i"%(cluster_id%10), alpha=0.4
        )
        if highlight_feature_residues:
            # Find the binding pocket (feature) residues
            pair = np.array([n.split('-') for n in features], dtype=int)
            bpid = np.unique(pair.flatten())
            isbp = [(r in bpid) for r in df_ca['resnum'] ]
            # Plot the RMSF of the binding pocket residues
            ax[cluster_id].bar(
                df_ca['resnum'][isbp], df_ca['RMSF'][isbp], 
                label=r'binding pocket C$\mathrm{\alpha}$', 
                width=1.0, color="C%i"%(cluster_id%10), alpha=1.0
            )
        # Legend
        ax[cluster_id].legend(framealpha=1)
    # Format and Labels
    for axi in ax:
        axi.set_xlim(np.min(resnum_min)-0.5, np.max(resnum_max)+0.5)
        axi.set_ylim([0,4.5])
        axi.set_ylabel('RMSF [$\mathrm{\AA}$]')
    ax[-1].set_xlabel('residue number')
    fig.tight_layout()
    # Save the figure as PNG and PDF files
    fig.savefig(out_plot+'.png', dpi=300)
    fig.savefig(out_plot+'.pdf', dpi=300)
    plt.close(fig)


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
    parser.add_argument('--ref_sel_align', dest='ref_asl_align', type=str, help='alignment selection for the reference file')
    parser.add_argument('--ref_sel_write', dest='ref_asl_write', type=str, help='output selection for the reference file')
    parser.add_argument('--align_avg', dest='align_avg', action='store_true', help='align the final average structure to the reference, using all atoms in the output selection.', default=False) 
    parser.add_argument('--threshold', dest='threshold', type=float, help='threshold RMSD to exclude frames', default=None)
    args = parser.parse_args()

    # Load selection files
    with open(args.sel_file) as sf:
        asl_align = [line.rstrip() for line in sf]
    with open(args.write_sel_file) as sf:
        asl_write = [line.rstrip() for line in sf]

    rmsf_per_atom, pos_average, cms_model_ref_new = calculate_rmsf(
        args.reference_fn, args.cms_files, args.trj_files, 
        args.ref_asl_align, args.ref_asl_write, asl_align, asl_write, 
        align_avg=args.align_avg, threshold=args.threshold)
   
    # Write the RMSF to a CSV file.  
    output = pd.DataFrame()
    output['RMSF'] = rmsf_per_atom
    output['pdbres'] = [a.pdbres for a in cms_model_ref_new.atom]
    output['resnum'] = [a.resnum for a in cms_model_ref_new.atom]
    output['pdbname'] = [a.pdbname for a in cms_model_ref_new.atom]
    output.to_csv(args.output_filename+'.csv')

    # Write the RMSF on reference structure.
    out_fn_ref = args.output_filename+'_rmsf_ref.cms'
    _ = write_coordinates(out_fn_ref, cms_model_ref_new, xyz=None, sigma=rmsf_per_atom)

    # Write the RMSF on average structure.
    out_fn_avg = args.output_filename+'_rmsf_avg.cms'
    _ = write_coordinates(out_fn_avg, cms_model_ref_new, pos_average, sigma=rmsf_per_atom)
    
