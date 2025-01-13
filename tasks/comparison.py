import os
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
import scipy.spatial
import scipy.spatial.distance
import matplotlib.pyplot as plt
from .io_features import read_features_from_csv_files, sort_features, get_feature_data


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


def plot_most_different_distributions(jsd_sorted, feat_i, feat_a, data_i, data_a, out_plot, showstart=False, feature_type='ca-distance'):

    fig, ax = plt.subplots(5,4, figsize=[10,8], dpi=300)
    ax = ax.flatten()

    for axi, fname in zip(ax, jsd_sorted[:20,0]):
        
        fdata_i = get_feature_data(feat_i, data_i, fname)
        fdata_a = get_feature_data(feat_a, data_a, fname)
        all_avg = np.mean(np.concatenate([fdata_i, fdata_a]))
        fdata_i = fdata_i - all_avg
        fdata_a = fdata_a - all_avg
        
        # Starting values
        if showstart:
            axi.axvline(fdata_i[0], lw=2, color='C0')
            axi.axvline(fdata_a[0], lw=2, color='C1')

        hist_c, bins_c = np.histogram(np.concatenate([fdata_i,fdata_a]), bins=60, density=True)
        hist_i, bins_i = np.histogram(fdata_i, bins=bins_c, density=True)
        hist_a, bins_a = np.histogram(fdata_a, bins=bins_c, density=True)
        bin_centers_c = .5 * (bins_c[1:] + bins_c[:-1])
        
        axi.plot(bin_centers_c, hist_i, lw=2)
        axi.plot(bin_centers_c, hist_a, lw=2)
        axi.fill_between(bin_centers_c, hist_i, alpha=0.5)
        axi.fill_between(bin_centers_c, hist_a, alpha=0.5)
        axi.fill_between(bin_centers_c, np.min([hist_i,hist_a], axis=0), alpha=0.25, color='k')
        
        axi.set_yticks([])
        axi.set_ylim(bottom=0)
        if feature_type == 'ca-distance':
            xlabel = r'distance CA %s - CA %s' % (fname.split('-')[0], fname.split('-')[1])
        else:
            xlabel = r'%s %s' % (feature_type, fname)
        axi.set_xlabel(xlabel)
        
    fig.tight_layout()  
    fig.savefig(out_plot+'.png', dpi=300)
    fig.savefig(out_plot+'.pdf', dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    """ 
    Perform systematic comparison of C-alpha distances.
    """
        
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='info_file', type=str, help='info csv file')
    parser.add_argument('-f', dest='feature_files', nargs='+', type=str, help='feature csv files')
    parser.add_argument('-s', dest='select_file', type=str, help='Name of the input file in csv format determining residue selections for binding pocket and restraints.') 
    parser.add_argument('-o', dest='output_dir', type=str, help='name of the output directory', default='comparison')   
    args = parser.parse_args()

    # Read the input files
    print("\n* - Reading the input files. - *\n")
    simulations = pd.read_csv(args.info_file)
    selections = pd.read_csv(args.select_file)
    names, data, origin, orig_id = read_features_from_csv_files(args.feature_files)

    # Split data in active and inactive (accoring to their origin simulation)
    act_origin = list(simulations[simulations['Start_Label']=='active'].index)
    ina_origin = list(simulations[simulations['Start_Label']=='inactive'].index)
    print('Inactive Simulations:', ina_origin, '\nActive Simulations:  ', act_origin, '\n')
    print('Shape of the total data:   ', data.shape)
    data_i = data[[o in ina_origin for o in origin]]
    data_a = data[[o in act_origin for o in origin]]
    print('Shape of the inactive data:', data_i.shape)
    print('Shape of the active data:  ', data_a.shape) 

    # Run the relative-entropy analysis
    data_names, jsd, kld_ab, kld_ba = relative_entropy_analysis(
        names, names, data_i, data_a, 
        bin_width=None, bin_num=10, verbose=True
        )
    # Sort the features by how much their distributions differ
    jsd_sorted = sort_features(data_names, jsd)
    out_data = pd.DataFrame(jsd_sorted, columns=['Distance','JSD'])
    out_csv = os.path.join(args.output_dir, 'ca-distances_sorted-by-jsd.csv')
    out_data.to_csv(out_csv)
    # Plot the 20 most different distributions
    out_file = os.path.join(args.output_dir, 'ca-distances_largest-jsd')
    plot_most_different_distributions(jsd_sorted, data_i, data_a, out_file)
    