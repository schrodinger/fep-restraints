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
    return msys_model, cms_model, trj


def copy_topology(reference, output_name):
    _, cms_model = topo.read_cms(reference)
    model = cms_model.copy()
    with structure.StructureWriter(output_name) as writer:
        writer.append(model.fsys_ct)
        for st in model.comp_ct:
            writer.append(st)


def write_frames(cms_model, trajectory, frame_numbers, out_dir, frame_names=None):
    
    # Determine the names of the frames (without extension)
    if frame_names is None:
        frame_names = { n:'frame_%04i'%n for n in frame_numbers }
    else:
        assert len(frame_names) == len(frame_numbers)
        frame_names = { n: frame_names[i] for i,n in enumerate(frame_numbers)}
        
    # Work with a copy of the model
    model = cms_model.copy()
    # Iterate through the trajectory
    for f, frame in enumerate(trajectory):
        # Only proceed if the current frame is among the selected
        if f in frame_numbers:
            # Update the ct positions to the traj positions.
            # We have to do this because traj includes pseudo-atoms.
            # And we have to do it for each ct separately to preserve information about which atom is in which ct.
            for i, ct in enumerate(model.comp_ct):
                ct_aids = model.get_fullsys_ct_atom_index_range(i)
                ct_gids = model.convert_to_gids(ct_aids, with_pseudo=False)
                topo.update_ct(ct, model, frame, is_fullsystem=False, allaid_gids=ct_gids)
            model.synchronize_fsys_ct()
            # Write the structure to the designated name
            out_fname = os.path.join(out_dir, frame_names[f])
            with structure.StructureWriter(out_fname+'.cms') as writer:
                writer.append(model.fsys_ct)
                for st in model.comp_ct:
                    writer.append(st)
    return


def extract_frames_by_value(trj_files, output_name, csv_files, value, property='Cluster_ID'):
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
    traj.write_traj(frame_list, output_name)
    