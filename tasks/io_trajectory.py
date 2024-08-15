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


def select_protein_residue(model, resnum, chain=None):
    """
    Selects a protein residue based on the residue number and optional chain.

    Parameters
    ----------
    model : object
        The molecular system model.
    resnum : int
        The residue number.
    chain : str, optional
        The chain identifier. Defaults to None.

    Returns
    -------
    object
        The selected residue object.

    """
    selstring = f"res. {resnum} AND at. CA"
    if chain is not None:
        selstring += f" AND chain. {chain}"
    ca = model.select_atom(selstring)
    if len(ca) == 0:
        print(f"Warning: No CA atom found for residue number {resnum}.")
        return None
    elif len(ca) > 1:
        print(f"Warning: More than one CA atom found for residue number {resnum}. Using the first one.")
    return model.atom[ca[0]].residue 
 

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


def extract_subset_model(cms_model, aid):
    num_atoms_all_cts = 0
    cms_model_new = cms_model.copy()
    for ct in cms_model_new.comp_ct:
        # Construct the index list of atoms to keep for this ct
        aidlist_write_ct = [_a-num_atoms_all_cts for _a in aid[num_atoms_all_cts:num_atoms_all_cts+ct.atom_total+1]]
        # Update total number of atoms across components
        num_atoms_all_cts += ct.atom_total
        # Construct the index list of atoms to delete for this ct
        aidlist_all = range(1, ct.atom_total+1)
        aidlist_del = list(set(aidlist_all) - set(aidlist_write_ct))
        # Delete the atoms from this ct
        ct.deleteAtoms(aidlist_del, renumber_map=True)
        # Delete all virtual and pseudo atoms from this ct
        pseudo_indices = [at.index for at in ct.ffio.pseudo]
        virtual_indices = [at.index for at in ct.ffio.virtual]
        ct.ffio.deletePseudos(pseudo_indices)
        ct.ffio.deleteVirtuals(virtual_indices)
    cms_model_new.synchronize_fsys_ct()
    return cms_model_new


def write_values_to_temperature_factor(cms_in, values, resnums, cms_out=None, chains=None, default_value=0.0):
    """
    Write values to the temperature factor of the residues in a structure, optionally write them to a CMS file.

    Parameters
    ----------
    cms_in : str
        The input CMS file. 
    values : list
        The values to be written.
    resnums : list
        The residue numbers.
    cms_out : str, optional
        The output CMS file. Defaults to None.
    chains : list, optional
        The chain identifiers. Defaults to None.
    default_value : float, optional
        The default value for the temperature factor. Defaults to 0.0.

    Returns
    -------
    object
        The structure object

    """
    st = structure.StructureReader.read(cms_in)
    if chains is None:
        chains = ['A'] * len(resnums)
    for res in st.residue:
        res.temperature_factor = default_value
    for v, r, c in zip(values, resnums, chains):
        res = st.findResidue(f"{c}:{r}")
        if res is not None:
            res.temperature_factor = v
    if cms_out:
        st.write(cms_out)
    return st


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


def extract_frames_by_value(top_files, trj_files, output_name, csv_files, value, property='Cluster_ID', 
                            start_frame=0, end_frame=None, step=1, asl_strings=None):
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
    out_cms_fname = output_name+'-out.cms'
    out_trj_fname = output_name+'.xtc'
    # Use all atoms if no asl strings are provided
    if asl_strings is None:
        asl_strings = ['all' for _ in range(len(trj_files))]
    # Check if the number of files is the same for all inputs
    assert len(trj_files) == len(csv_files) == len(top_files) == len(asl_strings)
    # Iterate through the files
    frame_list = []
    model_list = []
    for csv, top, trj, asl in zip(csv_files, top_files, trj_files, asl_strings):
        num_df = pd.read_csv(csv)
        clust_id = num_df[property]
        # Load the trajectory and CMS model
        msys_model, cms_model = topo.read_cms(top)
        print('Number of atoms in the system:', cms_model.atom_total)
        trajectory = traj.read_traj(trj)[start_frame:end_frame:step]
        print('Length of CSV file:', len(clust_id), ' Length of trajectory:', len(trajectory))
        # Select the atoms for output
        print('Output Selection:', str(asl))
        aidlist_write = cms_model.select_atom(str(asl)) 
        gidlist_write = topo.aids2gids(cms_model, aidlist_write, include_pseudoatoms=False)
        print("Selected %i atoms for output."%(len(gidlist_write)))
        # Construct a subsystem with the selected atoms
        out_cms_model = extract_subset_model(cms_model, aidlist_write)
        print('Number of atoms in the output:', out_cms_model.atom_total)
        # Extract frames with the given value
        assert len(clust_id) == len(trajectory)
        for clid, frame in zip(clust_id, trajectory):
            if clid == value:
                reduced_frame = frame.reduce(gidlist_write)
                frame_list.append(reduced_frame)
    # Write the extracted frames to a new trajectory
    out_cms_model.fix_filenames(out_cms_fname, out_trj_fname)
    out_cms_model.write(out_cms_fname)
    traj.write_traj(frame_list, out_trj_fname)
    