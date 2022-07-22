import os
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
import scipy.spatial
import scipy.spatial.distance
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from torch import _sobol_engine_initialize_state_
from utils.io_features import read_features_from_csv_files


def relative_entropy_analysis(features_a, features_b, all_data_a, all_data_b, bin_width=None, bin_num=10, verbose=True, override_name_check=False):
    """
    Calculates the Jensen-Shannon distance and the Kullback-Leibler divergences for each feature from two ensembles.
    
    Parameters
    ----------
        features_a : list of str
            Feature names of the first ensemble. 
            Can be obtained from features object via .describe().
        features_b : list of str
            Feature names of the first ensemble. 
            Can be obtained from features object via .describe().
            Must be the same as features_a. Provided as a sanity check. 
        all_data_a : float array
            Trajectory data from the first ensemble. Format: [frames,frame_data].
        all_data_b : float array
            Trajectory data from the second ensemble. Format: [frames,frame_data].
        bin_width : float, default=None
            Bin width for the axis to compare the distributions on. 
            If bin_width is None, bin_num (see below) bins are used and the width is determined from the common histogram.
        bin_num : int, default=10
            Number of bins for the axis to compare the distributions on (only if bin_width=None).
        verbose : bool, default=True
            Print intermediate results.
        override_name_check : bool, default=False
            Only check number of features, not their names.
    
    Returns
    -------
        data_names : list of str
            Feature names.
        data_jsdist : float array
            Jensen-Shannon distance for each feature.
        data_kld_ab : float array
            Kullback-Leibler divergences of data_a wrt to data_b.
        data_kld_ba : float array
            Kullback-Leibler divergences of data_b wrt to data_a.
        
    """
    all_data_a, all_data_b = all_data_a.T, all_data_b.T
    # Assert that the features are the same and data sets have same number of features
    if override_name_check:
        assert len(features_a) == len(features_b)
    else:
        assert features_a == features_b
    assert all_data_a.shape[0] == all_data_b.shape[0] 
    # Extract the names of the features
    data_names = features_a
    # Initialize relative entropy and average value
    data_jsdist = np.zeros(len(data_names))
    data_kld_ab = np.zeros(len(data_names))
    data_kld_ba = np.zeros(len(data_names))
    data_avg    = np.zeros(len(data_names))
    # Loop over all features
    for i in range(len(all_data_a)):
        data_a = all_data_a[i]
        data_b = all_data_b[i]
        # Combine both data sets
        data_both = np.concatenate((data_a,data_b))
        data_avg[i] = np.mean(data_both)
        # Get bin values for all histograms from the combined data set
        if bin_width is None:
            bins = bin_num
        else:
            bins_min = np.min( data_both )
            bins_max = np.max( data_both )
            bins = np.arange(bins_min,bins_max,bin_width)
        # Calculate histograms for combined and single data sets
        histo_both = np.histogram(data_both, bins = bins, density = True)
        histo_a = np.histogram(data_a, density = True, bins = histo_both[1])
        distr_a = histo_a[0] / np.sum(histo_a[0])
        histo_b = np.histogram(data_b, density = True, bins = histo_both[1])
        distr_b = histo_b[0] / np.sum(histo_b[0])
        # Calculate relative entropies between the two data sets (Kullback-Leibler divergence)
        data_kld_ab[i] = np.sum( sp.special.kl_div(distr_a,distr_b) )
        data_kld_ba[i] = np.sum( sp.special.kl_div(distr_b,distr_a) )
        # Calculate the Jensen-Shannon distance
        data_jsdist[i] = scipy.spatial.distance.jensenshannon(distr_a, distr_b, base=2.0)
        # Print information
        if verbose:
            print(i,'/',len(all_data_a),':', data_names[i]," %1.2f"%data_avg[i],
                  " %1.2f %1.2f %1.2f"%(data_jsdist[i],data_kld_ab[i],data_kld_ba[i]))   
    return data_names, data_jsdist, data_kld_ab, data_kld_ba


if __name__ == "__main__":
    """ 
    Perform systematic comparison of C-alpha distances.
    """
        
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_files', nargs='+', type=str, help='input csv files')
    parser.add_argument('-o', dest='output_base', type=str, help='base name of the output files', default='pca-clusters')
    parser.add_argument('-n', dest='n_components', type=int, help='number of top principal components to consider', default=3) 
    parser.add_argument('-k', nargs='+', dest='n_clusters', type=int, help='numbers k of clusters to attempt (arbitrary number)', default=[1,2,3,4])     
    parser.add_argument('-s', dest='random_state', type=int, help='random state for k-means algorithm', default=42)      
    parser.add_argument('-w', dest='write_pc', action='store_true', help='write coordinates in PC space', default=False)
    parser.add_argument('--out_names', dest='out_names', nargs='+', type=str, default=None) 
    args = parser.parse_args()

    names, data, origin, orig_id = read_features_from_csv_files(args.input_files)
    print(origin)

    # Print some information
    print("Origin:")
    for o in np.sort(np.unique(origin)):
        n_frames = np.sum(origin==o)
        top_file = args.input_files[o]
        print(" - %s: %5i frames"%(top_file, n_frames))

    