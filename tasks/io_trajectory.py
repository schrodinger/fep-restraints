from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import pandas as pd
import os.path


def get_an_aid(cms_model, asl_string, none_for_zero=False):
    """
    Get an atom identifier (aid) based on the given CMS model and ASL string.

    Parameters:
        cms_model (CMSModel): The CMS model object.
        asl_string (str): The ASL string used to select the atom.
        none_for_zero (bool, optional): If True, return None if no atom is found. 
            If False, raise an exception. Default is False.

    Returns:
        int or None: The atom identifier (aid) if found, or None if none_for_zero is True and no atom is found.

    Raises:
        Exception: If the number of chosen atoms is not equal to 1.
    """
    aid_list = cms_model.select_atom(asl_string)
    if (len(aid_list) == 0) and none_for_zero:
        return None
    elif (len(aid_list) != 1):
        raise Exception("The number of chosen atoms is %d."%(len(aid_list)))
    else:
        return(aid_list[0])
    

def load_trajectory(cms_file, xtc_file, trj_dir):
    """
    Load a molecular dynamics trajectory.

    Parameters
    ----------
    cms_file : str
        Path to the CMS file.
    xtc_file : str
        Path to the XTC file (optional).
    trj_dir : str
        Path to the trajectory directory (optional).

    Returns
    -------
    msys_model : object
        The molecular system model.
    cms_model : CMSModel 
        The CMS model.
    trj : object
        The loaded trajectory.

    """
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
    """
    Copy the topology from a reference file and save it to a new file.

    Parameters
    ----------
    reference : str
        The path to the reference file.
    output_name : str
        The name of the output file.

    """
    _, cms_model = topo.read_cms(reference)
    model = cms_model.copy()
    with structure.StructureWriter(output_name) as writer:
        writer.append(model.fsys_ct)
        for st in model.comp_ct:
            writer.append(st)


def write_frames(cms_model, trajectory, frame_numbers, out_dir, frame_names=None):
    """
    Write selected frames from a trajectory to separate structure files.

    Parameters
    ----------
    cms_model : CmsModel
        The CMS model object.
    trajectory : list
        The trajectory containing the frames.
    frame_numbers : list
        The indices of the frames to be written.
    out_dir : str
        The output directory where the structure files will be saved.
    frame_names : dict, optional
        A dictionary mapping frame indices to custom frame names. 
        If not provided, default frame names will be used.

    Returns
    -------
    None

    """
    
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


def extract_frames_by_value(trj_files, output_name, csv_files, value, property='Cluster_ID', start_frame=0, end_frame=None, step=1):
    """
    Extracts frames from trajectory files based on a given value of a property.

    Parameters:
    -----------
        trj_files : list
            List of trajectory file paths.
        output_name : str
            Output file name for the extracted frames.
        csv_files : list
            List of CSV file paths.
        value : int or str
            Value of the property to filter frames by.
        property : str, optional
            Name of the property to filter frames by. Defaults to 'Cluster_ID'.
        start_frame : int, optional
            Index of the first frame from the trajectories. Defaults to 0.
        end_frame : int, optional
            Index of the last frame from the trajectories. Defaults to None, which extracts all frames.
        step : int, optional
            Step size for trajectories. Defaults to 1.

    Returns:
    --------
    None

    """
    frame_list = []
    for csv, trj in zip(csv_files, trj_files):
        num_df = pd.read_csv(csv)
        clust_id = num_df[property]
        trajectory = traj.read_traj(trj)[start_frame:end_frame:step]
        print('Length of CSV file:', len(clust_id), ' Length of trajectory:', len(trajectory))
        assert len(clust_id) == len(trajectory)
        for clid, frame in zip(clust_id, trajectory):
            if clid == value:
                frame_list.append(frame)
    traj.write_traj(frame_list, output_name)
    