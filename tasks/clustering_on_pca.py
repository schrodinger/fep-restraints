import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from torch import _sobol_engine_initialize_state_
from utils.io_features import read_features_from_csv_files


def kmeans_on_pca(pc, k, rs, origin, orig_id, output_base, input_files=None, write_pc=True, out_names=None):

    # Define the algorithm
    kmeans = KMeans(n_clusters = k, init = 'k-means++', random_state = rs)
    # Fit and transform the data to cluster-distance space.
    cdist = kmeans.fit_transform(pc.T)

    # Determine clusters and their centers.
    clusters = kmeans.fit_predict(pc.T)
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
        if write_pc:
            for i, pci in enumerate(pc):
                output['PC%02i'%(i+1)] = pci[is_from_origin]
        output.to_csv(cl_file_name)
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
    return cids, sizes, cc_orig_sim, cc_orig_id, kmeans.inertia_, cl_files, summary_name 



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
    print(origin)

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
        cids, sizes, cc_orig_sim, cc_orig_id, inertia, cl_files, sum_file = kmeans_on_pca(
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
