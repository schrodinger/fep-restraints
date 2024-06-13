import numpy as np
import pandas as pd
from schrodinger.application.desmond.packages import analysis


def write_features_to_csv(out_file, features_data, feature_names, frame_time):
    """
    Write features to a CSV file.

    Parameters
    ----------
    out_file : str
        File name for the output CSV file.
    feature_names : list of str
        Names of the features.
    features_data : numpy array
        Data for the features. Format: [features, frames].
    frame_time : numpy array
        Time of each frame in ps.
    """ 
    # create and write csv file
    results = pd.DataFrame()
    results['time (ps)'] = frame_time
    if (len(feature_names) > 1):  # 'features' is a two dimensional array
        for name, values in zip(feature_names, features_data):
            results[name] = values
    else:
        results[feature_names[0]] = features_data
    results.to_csv(out_file, index=False)


def read_features_from_csv(csv_file):
    """
    Load features from a CSV file.

    Parameters
    ----------
    csv_file : str
        File name for the input CSV file.

    Returns
    -------
    feature_names : list of str
        Names of the features.
    features_data : numpy array
        Data for the features. Format: [frames, frame_data].

    """ 
    df = pd.read_csv(csv_file) 
    feature_names = list(df.keys()[1:]) # start from 1 to skip the time column
    feature_data = np.zeros([len(df),len(feature_names)])
    for i,f in enumerate(feature_names):
        feature_data[:,i] = df[f]
    return feature_names, feature_data 


def read_features_from_csv_files(csv_names):
    """
    Load features from multiple CSV files.

    Parameters
    ----------
    csv_files : list of str
        File names for the input CSV files.

    Returns
    -------
    feature_names : list of str
        Names of the features.
    features_data : numpy array
        Data for the features. Format: [frames, frame_data].

    """ 
    all_data = []
    orig_sim = []
    orig_idx = []
    old_feat = None
    for i,csv in enumerate(csv_names):
        new_feat, new_data = read_features_from_csv(csv)
        all_data.append(new_data)
        orig_sim += [i for j in range(len(new_data))]
        orig_idx += [j for j in range(len(new_data))]
        if i > 0 and old_feat != new_feat:
            print('Warning: Inconsistent feature names!')
        old_feat = new_feat.copy()
    all_data = np.concatenate(all_data)
    orig_sim = np.array(orig_sim, dtype=int)
    orig_idx = np.array(orig_idx, dtype=int)
    return new_feat, all_data, orig_sim, orig_idx


def sort_features(names, sortby):
    """
    Sorts features by a list of values.
    
    Parameters
    ----------
    names : str array
        Array of feature names.
    sortby : float array
        Array of the values to sort the names by.
        
    Returns
    -------
    sort : array of tuples [str,float]
        Array of sorted tuples with feature and value.
    
    """
    # Get the indices of the sorted order
    sort_id = np.argsort(sortby)[::-1]  
    # Bring the names and values in the right order
    sorted_names  = []
    sorted_values = []
    for i in sort_id:
        sorted_names.append(np.array(names)[i])
        sorted_values.append(sortby[i])
    sn, sv = np.array(sorted_names), np.array(sorted_values)
    # Format for output
    sort = np.array([sn,sv]).T
    return sort


def get_feature_data(feat, data, feature_name):
    """
    Returns the timeseries of one particular feature.
    
    Parameters
    ----------
        feat : list of str
            List with all feature names.
        data : float array
            Feature values data from the simulation.
        feature_name : str
            Name of the selected feature.
    
    Returns
    -------
        timeseries : float array
            Value of the feature for each frame.
    
    """
    # Select the feature and get its index.
    index = np.where( np.array( feat ) == feature_name )[0][0]
    # Extract the timeseries.
    timeseries = data[:,index]
    return timeseries
     

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


def calculate_ca_distances(msys_model, cms_model, tr, chain_ids, residue_numbers, residue_names=None, chain_id_in_name=False):
    """
    Calculates distances between pairs of C-alpha atoms from a trajectory.
    
    Parameters
    ----------
        msys_model : msys
            System model object from the simulation's topology.
        cms_model : cms
            CMS model object from the simulation's topology.
        tr : trajectory
            Trajectory object from the simulation.
        chain_ids : list
            List with the names of all the chains with selected atoms.
        residue_numbers : list
            List with lists of residue numbers for each chain.
        residue_names : list, optional
            List with list of residue names for each chain.

    Returns
    -------
        frame_time : float array
            Time for each frame.
        distance_names : list
            Names of all the distances.
        distances : list
            Data for the distances. Format: [features, frames].
    
    """
    assert len(chain_ids) == len(residue_numbers)
    # time
    frame_time = []
    for item in tr:
        frame_time.append(item.time)
    # define residue names
    if residue_names is None:
        residue_names = residue_numbers
    else:
        assert len(residue_names) == len(residue_numbers)
        for c, chain_id in enumerate(chain_ids):
            assert len(residue_names[c]) == len(residue_numbers[c])
    # define residues
    residues = []
    for c, chain_id in enumerate(chain_ids):
        for number, name in zip(residue_numbers[c], residue_names[c]):
            new_res = {'chain':chain_id, 'number':number, 'name':name}
            residues.append(new_res)
    # define analyzers
    analyzers = []
    distance_names = []
    for i in range(len(residues)):
        for j in range(i+1, len(residues)):
            chain_i, rnum_i, rname_i  = residues[i]['chain'], residues[i]['number'], residues[i]['name']
            chain_j, rnum_j, rname_j  = residues[j]['chain'], residues[j]['number'], residues[j]['name']
            first_asl = '(res.num %i) AND (atom.ptype " CA ") AND (chain.name %s)'%(rnum_i, chain_i)
            first_aid = get_an_aid(cms_model, first_asl)
            second_asl = '(res.num %i) AND (atom.ptype " CA ") AND (chain.name %s)'%(rnum_j, chain_j)
            second_aid = get_an_aid(cms_model, second_asl)
            analyzers.append(analysis.Distance(msys_model, cms_model, first_aid, second_aid))
            if chain_id_in_name:
                distance_names.append("%s:%s-%s:%s"%(chain_i, rname_i, chain_j, rname_j))
            else:
                distance_names.append("%s-%s"%(rname_i, rname_j))
    #compute result
    distances = analysis.analyze(tr, *analyzers)
    return frame_time, distance_names, distances


# BACKBONE TORSION DEFINITIONS
bb_torsion_ptypes = {}
bb_torsion_relres = {}
#  PHI (φ):   C(i-1)-N(i)-CA(i)-C(i)
bb_torsion_ptypes['PHI'] = ["C", "N", "CA", "C"]
bb_torsion_relres['PHI'] = [-1, 0, 0, 0]
#  PSI (ψ):   N(i)-CA(i)-C(i)-N(i + 1)
bb_torsion_ptypes['PSI'] = ["N", "CA", "C", "N"]
bb_torsion_relres['PSI'] = [0, 0, 0, 1]
#  OMEGA (ω): CA(i)-C(i)-N(i + 1)-CA(i + 1)
bb_torsion_ptypes['OMEGA'] = ["CA", "C", "N", "CA"]
bb_torsion_relres['OMEGA'] = [0, 0, 1, 1]


def calculate_backbone_torsions(msys_model, cms_model, tr, chain_ids, residue_numbers, residue_names=None, chain_id_in_name=False):
    """
    Calculates protein backbone torsions from a trajectory.
    
    Parameters
    ----------
        msys_model : msys
            System model object from the simulation's topology.
        cms_model : cms
            CMS model object from the simulation's topology.
        tr : trajectory
            Trajectory object from the simulation.
        chain_ids : list
            List with the names of all the chains with selected atoms.
        residue_numbers : list
            List with lists of residue numbers for each chain.
        residue_names : list, optional
            List with list of residue names for each chain.

    Returns
    -------
        frame_time : float array
            Time for each frame.
        torsion_names : list
            Names of all the torsions.
        torsions : list
            Data for the torsions. Format: [features, frames].
    
    """
    assert len(chain_ids) == len(residue_numbers)
    # Read the time for each frame
    frame_time = []
    for item in tr:
        frame_time.append(item.time)
    # Define the residue names
    if residue_names is None:
        residue_names = residue_numbers
    else:
        assert len(residue_names) == len(residue_numbers)
        for c, chain_id in enumerate(chain_ids):
            assert len(residue_names[c]) == len(residue_numbers[c])
    # Define the residues
    residues = []
    for c, chain_id in enumerate(chain_ids):
        for number, name in zip(residue_numbers[c], residue_names[c]):
            new_res = {'chain':chain_id, 'number':number, 'name':name}
            residues.append(new_res)
    # Define the analyzers
    analyzers = []
    torsion_names = []
    # In each backbone torsion, all atoms have the same chain but some have a different residue number/name
    for r in range(len(residues)):
        chain, rnum, rname  = residues[r]['chain'], residues[r]['number'], residues[r]['name']
        # There are multiple possible sidechain torsions, try them all!
        # Loop over each of the possible types (chi1-5)
        for torsion_type in bb_torsion_ptypes.keys():
            # Within each type, get the ptypes (atom names) and relative residue numbers
            atom_ptypes = bb_torsion_ptypes[torsion_type]
            atom_relres = bb_torsion_relres[torsion_type]
            ptype_i, ptype_j, ptype_k, ptype_l = atom_ptypes
            relrn_i, relrn_j, relrn_k, relrn_l = atom_relres
            # Construct the selection strings
            asl_i = '(res.num %i) AND (atom.ptype "%s") AND (chain.name %s)'%(rnum + relrn_i, ptype_i, chain)
            asl_j = '(res.num %i) AND (atom.ptype "%s") AND (chain.name %s)'%(rnum + relrn_j, ptype_j, chain)
            asl_k = '(res.num %i) AND (atom.ptype "%s") AND (chain.name %s)'%(rnum + relrn_k, ptype_k, chain)
            asl_l = '(res.num %i) AND (atom.ptype "%s") AND (chain.name %s)'%(rnum + relrn_l, ptype_l, chain)
            # Get the corresponding atom IDs (None if atom not in system)
            aid_i = get_an_aid(cms_model, asl_i, none_for_zero=True)
            aid_j = get_an_aid(cms_model, asl_j, none_for_zero=True)
            aid_k = get_an_aid(cms_model, asl_k, none_for_zero=True)
            aid_l = get_an_aid(cms_model, asl_l, none_for_zero=True)
            # Only add a torsion if all atoms are in the structure
            if aid_i and aid_j and aid_k and aid_l:
                analyzers.append(analysis.Torsion(msys_model, cms_model, aid_i, aid_j, aid_k, aid_l))
                if chain_id_in_name:
                    torsion_names.append("%s:%s-%s"%(chain, rname, torsion_type))
                else:
                    torsion_names.append("%s-%s"%(rname, torsion_type))
    # Compute the result
    torsions = analysis.analyze(tr, *analyzers)
    return frame_time, torsion_names, torsions


# SIDECHAIN TORSION DEFINITIONS
sc_torsion_ptypes = {}
sc_torsion_ptypes['CHI1'] = [
    ["N", "CA", "CB", "CG"],
    ["N", "CA", "CB", "CG1"],
    ["N", "CA", "CB", "SG"],
    ["N", "CA", "CB", "OG"],
    ["N", "CA", "CB", "OG1"]
]
sc_torsion_ptypes['CHI2'] = [
    ["CA", "CB", "CG", "CD"],
    ["CA", "CB", "CG", "CD1"],
    ["CA", "CB", "CG1", "CD"],
    ["CA", "CB", "CG1", "CD1"],
    ["CA", "CB", "CG", "OD1"],
    ["CA", "CB", "CG", "ND1"],
    ["CA", "CB", "CG", "SD"]
]
sc_torsion_ptypes['CHI3'] = [
    ["CB", "CG", "CD", "NE"],
    ["CB", "CG", "CD", "CE"],
    ["CB", "CG", "CD", "OE1"],
    ["CB", "CG", "SD", "CE"]
]
sc_torsion_ptypes['CHI4'] = [
    ["CG", "CD", "NE", "CZ"],
    ["CG", "CD", "CE", "NZ"]
]
sc_torsion_ptypes['CHI5'] = [
    ["CD", "NE", "CZ", "NH1"]
]


def calculate_sidechain_torsions(msys_model, cms_model, tr, chain_ids, residue_numbers, residue_names=None, chain_id_in_name=False):
    """
    Calculates protein sidechain torsions from a trajectory.
    
    Parameters
    ----------
        msys_model : msys
            System model object from the simulation's topology.
        cms_model : cms
            CMS model object from the simulation's topology.
        tr : trajectory
            Trajectory object from the simulation.
        chain_ids : list
            List with the names of all the chains with selected atoms.
        residue_numbers : list
            List with lists of residue numbers for each chain.
        residue_names : list, optional
            List with list of residue names for each chain.

    Returns
    -------
        frame_time : float array
            Time for each frame.
        torsion_names : list
            Names of all the torsions.
        torsions : list
            Data for the torsions. Format: [features, frames].
    
    """
    assert len(chain_ids) == len(residue_numbers)
    # Read the time for each frame
    frame_time = []
    for item in tr:
        frame_time.append(item.time)
    # Define the residue names
    if residue_names is None:
        residue_names = residue_numbers
    else:
        assert len(residue_names) == len(residue_numbers)
        for c, chain_id in enumerate(chain_ids):
            assert len(residue_names[c]) == len(residue_numbers[c])
    # Define the residues
    residues = []
    for c, chain_id in enumerate(chain_ids):
        for number, name in zip(residue_numbers[c], residue_names[c]):
            new_res = {'chain':chain_id, 'number':number, 'name':name}
            residues.append(new_res)
    # Define the analyzers
    analyzers = []
    torsion_names = []
    # In each sidechain torsion, all atoms have the same chain and residue number/name
    for r in range(len(residues)):
        chain, rnum, rname  = residues[r]['chain'], residues[r]['number'], residues[r]['name']
        # There are multiple possible sidechain torsions, try them all!
        # Loop over each of the possible types (chi1-5)
        for torsion_type in sc_torsion_ptypes.keys():
            # Within each type, loop over the potential sets of ptypes (atom names)
            for atom_ptypes in sc_torsion_ptypes[torsion_type]:
                ptype_i, ptype_j, ptype_k, ptype_l = atom_ptypes
                # Construct the selection strings
                asl_i = '(res.num %i) AND (atom.ptype "%s") AND (chain.name %s)'%(rnum, ptype_i, chain)
                asl_j = '(res.num %i) AND (atom.ptype "%s") AND (chain.name %s)'%(rnum, ptype_j, chain)
                asl_k = '(res.num %i) AND (atom.ptype "%s") AND (chain.name %s)'%(rnum, ptype_k, chain)
                asl_l = '(res.num %i) AND (atom.ptype "%s") AND (chain.name %s)'%(rnum, ptype_l, chain)
                # Get the corresponding atom IDs (None if atom not in system)
                aid_i = get_an_aid(cms_model, asl_i, none_for_zero=True)
                aid_j = get_an_aid(cms_model, asl_j, none_for_zero=True)
                aid_k = get_an_aid(cms_model, asl_k, none_for_zero=True)
                aid_l = get_an_aid(cms_model, asl_l, none_for_zero=True)
                # Only add a torsion if all atoms are in the structure
                if aid_i and aid_j and aid_k and aid_l:
                    analyzers.append(analysis.Torsion(msys_model, cms_model, aid_i, aid_j, aid_k, aid_l))
                    if chain_id_in_name:
                        torsion_names.append("%s:%s-%s"%(chain, rname, torsion_type))
                    else: 
                        torsion_names.append("%s-%s"%(rname, torsion_type))
    # Compute the result
    torsions = analysis.analyze(tr, *analyzers)
    return frame_time, torsion_names, torsions
