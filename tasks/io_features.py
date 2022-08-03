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
     

def get_an_aid(cms_model, asl_string):
    aid_list = cms_model.select_atom(asl_string)
    if (len(aid_list) != 1):
        raise Exception("The number of chosen atoms is %d."%(len(aid_list)))
    else:
        return(aid_list[0])


def calculate_ca_distances(msys_model, cms_model, tr, chain_ids, residue_numbers, residue_names=None):
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
            distance_names.append("%s-%s"%(rname_i, rname_j))
    #compute result
    distances = analysis.analyze(tr, *analyzers)
    return frame_time, distance_names, distances
