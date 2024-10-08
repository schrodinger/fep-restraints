import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.stats import norm
from torch import _sobol_engine_initialize_state_
from .io_features import read_features_from_csv_files
from .io_trajectory import write_values_to_temperature_factor 


def kmeans_on_pca(data_pca, k, rs, origin, orig_id, output_base, input_files=None, write_pc=True, out_names=None):

    # Define the algorithm
    kmeans = KMeans(n_clusters = k, init = 'k-means++', random_state = rs)
    # Fit and transform the data to cluster-distance space.
    cdist = kmeans.fit_transform(data_pca)

    # Determine clusters and their centers.
    clusters = kmeans.fit_predict(data_pca)
    centers = kmeans.cluster_centers_
    sizes = [np.sum(clusters==ic) for ic in range(k)]
    cids = np.arange(k)

    # Find centroids and their indices and origin files.
    centroids = np.argmin(cdist, axis=0)
    cc_orig_sim = origin[centroids]
    cc_orig_id = orig_id[centroids]

    # Write information about these clusters and their components to a csv file
    cl_files = []
    for i, ii in enumerate(np.unique(origin)):
        if out_names is None:
            cl_file_name = output_base + '_sim%04i.csv'%ii
        else:
            cl_file_name = output_base + '_' + out_names[i] + '.csv'
        print('Writing cluster information to', cl_file_name) 
        is_from_origin = origin==ii
        output = pd.DataFrame()
        output['Frame_ID'] = np.arange(np.sum(is_from_origin))
        output['Cluster_ID'] = clusters[is_from_origin]
        for j in range(k):
            output['Distance_to_Cluster_%02i'%j] = cdist[is_from_origin,j]
        if write_pc:
            for j, pcj in enumerate(data_pca.T):
                output['PC%02i'%(j+1)] = pcj[is_from_origin]
        output.to_csv(cl_file_name, index=False)
        cl_files.append(cl_file_name)

    # Write summary information
    summary_name = output_base + '_summary.csv'
    print('Writing summary information to', summary_name)
    output = pd.DataFrame()
    results = [cids, sizes, cc_orig_sim, cc_orig_id]
    columns = ['Cluster_ID', 'Size', 'COrig_File_ID', 'Centroid_Original_Index']
    if write_pc:
        for i in range(len(centers[0])):
            columns += ['Center_PC%02i'%(i+1)]
            results += [centers[:,i]]
    for col, res in zip(columns, results):
        output[col] = res
    if input_files is not None:
        cc_orig_fn = [input_files[o] for o in cc_orig_sim]
        output['Centroid_Original_File'] = cc_orig_fn 
    output.to_csv(summary_name, index=False)
    print(output)
    
    # Return the sum of squared distances for this k.
    return cids, sizes, centers, cdist, cc_orig_sim, cc_orig_id, kmeans.inertia_, cl_files, summary_name 


def scatterplot_pca_by_system(data_pca, index1, index2, origin, simulations, out_file, showstart=False):
    """
    PC indices are zero-based!
    """
    systems = simulations['System_Name'].unique()
    pc1 = data_pca[:,index1]
    pc2 = data_pca[:,index2]
    fig, ax = plt.subplots(1,1, figsize=[3,3], dpi=300)
    means = []
    starts = []
    for sys_id, sys in enumerate(systems):
        # Find PC values of all data points from each system
        sys_origin = list(simulations[simulations['System_Name']==sys].index)
        sys_pc1 = pc1[[o in sys_origin for o in origin]]
        sys_pc2 = pc2[[o in sys_origin for o in origin]]
        # Plot all data points of this system
        ax.plot(sys_pc1, sys_pc2, '.', mew=0, ms=2, alpha=0.2, color='C%i'%sys_id)
        means.append([np.mean(sys_pc1), np.mean(sys_pc2)])
        starts.append([sys_pc1[0], sys_pc2[0]])
    for sys_id, sys in enumerate(systems):
        if showstart: 
            ax.plot(*starts[sys_id], 's', mew=1, mec='k', alpha=1, color='C%i'%sys_id) 
        ax.plot(*means[sys_id], 'o', mew=1, mec='k', alpha=1, color='C%i'%sys_id, label=sys)
    # Format and labels
    ax.set_xlim(np.min(data_pca[:,index1]), np.max(data_pca[:,index1]))
    ax.set_ylim(np.min(data_pca[:,index2]), np.max(data_pca[:,index2]))
    ax.set_xlabel('PC%i'%(index1+1))
    ax.set_ylabel('PC%i'%(index2+1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_file+'.pdf', dpi=300)
    fig.savefig(out_file+'.png', dpi=300)
    plt.close(fig)


def plot_pca_by_system(data_pca, origin, simulations, out_file, showstart=False):  
    systems = simulations['System_Name'].unique()
    for i, pci in enumerate(data_pca.T):
        # Calculate the general bins and their centers
        hist_c, bins_c = np.histogram(pci, bins=50)
        bin_centers = .5 * (bins_c[1:] + bins_c[:-1])
        # Start the figure
        fig, ax = plt.subplots(1, 1, figsize=[3,2], dpi=300)
        # Loop over systems
        for sys_id, sys in enumerate(systems):
            # Find PC values of all data points from each system
            sys_origin = list(simulations[simulations['System_Name']==sys].index)
            sys_pci = pci[[o in sys_origin for o in origin]]    
            if showstart:
                ax.axvline(sys_pci[0], lw=2, color='C%i'%sys_id)                
            # Calculate the histogram on the general bins and plot it
            hist_i, bins_i = np.histogram(sys_pci, bins=bins_c)
            ax.plot(bin_centers, hist_i, lw=1, alpha=1, label=sys)
            ax.fill_between(bin_centers, hist_i, alpha=0.25)
        # Format and labels
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('PC%i'%(i+1))
        ax.set_ylim(bottom=0)    
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_file+'_PC%02i.pdf'%(i+1), dpi=300)
        fig.savefig(out_file+'_PC%02i.png'%(i+1), dpi=300)
        plt.close(fig)


def elbow_plot(num_clusters, sum_squ_dist, out_file):
    fig,ax = plt.subplots(1, 1, figsize=[3,2.25], dpi=300)
    ax.plot(num_clusters, sum_squ_dist, 'o')
    ncmin, ncmax = np.min(num_clusters), np.max(num_clusters)
    ssmin, ssmax = np.min(sum_squ_dist), np.max(sum_squ_dist)
    # Format and labels
    ax.set_xticks(np.arange(ncmin,ncmax+1))
    ax.set_xticklabels(np.arange(ncmin,ncmax+1), size=7)
    ax.set_yticks([])
    ax.set_xlim([ncmin-0.5,9.5])
    ax.set_ylim([0,1.1*ssmax])
    ax.set_xlabel('number of clusters')
    ax.set_ylabel('sum of sq. dist.')
    fig.tight_layout()
    fig.savefig(out_file+'.png', dpi=300)
    fig.savefig(out_file+'.pdf', dpi=300)
    plt.close(fig)


def pc_cluster_plot(data_pca, cl_files_k, centers, out_pca_cl, index1=0, index2=1):
    pc1 = data_pca[:,index1]
    pc2 = data_pca[:,index2]
    # Get the data 
    cluster_data = []
    for clf in cl_files_k:
        new_data = pd.read_csv(clf)
        cluster_data.append(new_data['Cluster_ID'])
    cluster_data = pd.concat(cluster_data)
    k = len(np.unique(cluster_data))
    # Plot PCA with clusters
    fig,ax = plt.subplots(1, 1, figsize=[3,3], dpi=300)
    for cluster_id in range(k):
        # Find PC values of all data points from each system
        is_in_cluster = [c == cluster_id for c in cluster_data]
        cl_pc1 = pc1[is_in_cluster]
        cl_pc2 = pc2[is_in_cluster]
        ax.plot(cl_pc1, cl_pc2, '.', mew=0, ms=2, alpha=0.2, color="C%i"%(cluster_id%10))
    # Plot the cluster centers
    for cluster_id in range(k):
        cc1, cc2 = centers[cluster_id,index1], centers[cluster_id,index2]
        ax.plot(cc1, cc2, 'o', alpha=0.5, mec='k', mew=1,
            color="C%i"%(cluster_id%10), label='%i'%cluster_id)
    # Format and labels
    ax.set_xlim(np.min(data_pca[:,index1]), np.max(data_pca[:,index1]))
    ax.set_ylim(np.min(data_pca[:,index2]), np.max(data_pca[:,index2]))
    ax.set_xlabel('PC%i'%(index1+1))
    ax.set_ylabel('PC%i'%(index2+1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(fontsize=8) 
    fig.tight_layout()
    fig.savefig(out_pca_cl+'.pdf', dpi=300)
    fig.savefig(out_pca_cl+'.png', dpi=300)
    plt.close(fig)


def calculate_residue_contributions(csv, feature_type):
    """
    Calculate the maximum contribution of each residue to each principal component.
    
    Parameters
    ----------
    csv : str
        The path to the CSV file containing the features.
    feature_type : str
        The type of features to consider. Possible values are 'ca-distance', 'bb-torsion', or 'sc-torsion'.
    
    Returns
    -------
    residue_contributions : pandas.DataFrame
        A DataFrame containing the maximum contribution of each residue to each principal component.
        The 'Residue' column contains the unique residue numbers, and each principal component is represented by a separate column.
    """
    
    # Load the features
    df = pd.read_csv(csv)
    features = df['Feature']
    pc_keys = [key for key in df.keys() if key.startswith('PC')]
    
    # Construct the data frame for output
    residue_contributions = pd.DataFrame()
    if feature_type == 'ca-distance':
        pair = np.array([n.split('-') for n in features], dtype=int)
        df['Residue 1'] = pair[:,0]
        df['Residue 2'] = pair[:,1]
        unique_residue_numbers = np.unique(pair.flatten())
        residue_contributions['Residue'] = unique_residue_numbers
        # Calculate the maximum contribution of all the features in each residue
        for i, pc in enumerate(pc_keys):
            residue_contributions_pc = np.zeros(len(unique_residue_numbers))
            for j, num in enumerate(unique_residue_numbers):
                df_res1 = df[df['Residue 1']==num]
                df_res2 = df[df['Residue 2']==num]
                df_res = df_res2.append(df_res1)
                residue_contributions_pc[j] = np.max(np.abs(df_res[pc]))
            residue_contributions[pc] = residue_contributions_pc     
    elif feature_type == 'bb-torsion' or feature_type == 'sc-torsion':
        residue_numbers = np.array([n.split('-')[0] for n in features], dtype=int)
        df['Residue'] = residue_numbers
        unique_residue_numbers = np.unique(residue_numbers)
        residue_contributions['Residue'] = unique_residue_numbers
        # Calculate the maximum contribution of all the features in each residue
        for i, pc in enumerate(pc_keys):
            residue_contributions_pc = np.zeros(len(unique_residue_numbers))
            for j, num in enumerate(unique_residue_numbers):
                df_res = df[df['Residue']==num]
                residue_contributions_pc[j] = np.max(np.abs(df_res[pc]))
            residue_contributions[pc] = residue_contributions_pc 
            
    return residue_contributions


def plot_residue_contributions(residue_contributions, out_name):
    """
    Plot the residue contributions to principal components.
    """
    pc_keys = [key for key in residue_contributions.keys() if key.startswith('PC')]
    residues = residue_contributions['Residue']
    fig, ax = plt.subplots(len(pc_keys), 1, figsize=[6, 3 * len(pc_keys)], dpi=300)
    for i, pc in enumerate(pc_keys):
        ax[i].bar(residues, residue_contributions[pc])
        ax[i].set_xlim(np.min(residues), np.max(residues))
        ax[i].set_ylabel('max. abs. contrib. to PC%i' % (i + 1))
        fig.tight_layout()
    ax[-1].set_xlabel('residue number')
    fig.tight_layout()
    fig.savefig(out_name + '.pdf', dpi=300)
    fig.savefig(out_name + '.png', dpi=300)
    plt.close(fig)
    return


if __name__ == "__main__":
    """ 
    Perform k-means clustering in principal component space.
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

    # Print some information
    print("Origin:")
    for o in np.sort(np.unique(origin)):
        n_frames = np.sum(origin==o)
        top_file = args.input_files[o]
        print(" - %s: %5i frames"%(top_file, n_frames))

    # Perform PCA on the data.
    pca = PCA(n_components=args.n_components)
    pca.fit(data.T)
    pc = pca.components_
    print('PC shape:', pc.shape)
    print(pca.explained_variance_ratio_) 

    # Perform k-means clustering in PC space for various k values.
    sum_sqrd = []
    for ik, k in enumerate(args.n_clusters):  
        outputf = args.output_base+'_n%02i_s%02i_k%02i'%(args.n_components, args.random_state, k)
        cids, sizes, centers, cdist, cc_orig_sim, cc_orig_id, inertia, cl_files, sum_file = kmeans_on_pca(
            pc, k, args.random_state, origin, orig_id, output_base=outputf, 
            input_files=args.input_files, write_pc=args.write_pc, out_names=args.out_names
            )
        sum_sqrd.append(inertia)

    # Write information about all k values in this study
    file_name_ssd = args.output_base 
    file_name_ssd += '_n%02i_s%02i_ssd.csv'%(args.n_components, args.random_state)
    print('Writing the sums of the squared distances to', file_name_ssd)
    ssd = pd.DataFrame()
    ssd['Num_Clusters'] = args.n_clusters
    ssd['Sum_Squ_Dist'] = sum_sqrd 
    ssd.to_csv(file_name_ssd, index=False)
