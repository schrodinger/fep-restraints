import numpy as np
import pandas as pd


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
    df = pd.read_csv(csv_file) # skip the index column
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
    old_feat = None
    for i,csv in enumerate(csv_names):
        new_feat, new_data = read_features_from_csv(csv)
        all_data.append(new_data)
        if i > 0 and old_feat != new_feat:
            print('Warning: Inconsistent feature names!')
        old_feat = new_feat.copy()
    all_data = np.concatenate(all_data)
    return new_feat, all_data