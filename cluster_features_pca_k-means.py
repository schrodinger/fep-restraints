import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from torch import _sobol_engine_initialize_state_
from io_features import read_features_from_csv_files

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

        # Define the algorithm
        kmeans = KMeans(
            n_clusters = k, 
            init = 'k-means++', 
            random_state = args.random_state
        )
        # Fit and transform the data to cluster-distance space.
        cdist = kmeans.fit_transform(pc.T)
        # Calculate the sum of squared distances for this k.
        sum_sqrd.append(kmeans.inertia_) 

        # Determine clusters and their centers.
        clusters = kmeans.fit_predict(pc.T)
        centers = kmeans.cluster_centers_
        sizes = [np.sum(clusters==ic) for ic in range(k)]
        cids = np.arange(k)

        # Find centroids and their indices and origin files.
        centroids = np.argmin(cdist, axis=0)
        cc_orig_sim = origin[centroids]
        cc_orig_id = orig_id[centroids]
        cc_orig_fn = [args.input_files[o] for o in cc_orig_sim]

        # Write information about these clusters and their components to a csv file
        for ii, infile in enumerate(args.input_files):
            cl_file_name = args.output_base
            cl_file_name += '_n%02i_s%02i_k%02i_%s'%(args.n_components, args.random_state, k, os.path.basename(infile))
            print('Writing cluster information to', cl_file_name) 
            is_from_origin = origin==ii
            output = pd.DataFrame()
            output['Frame_ID'] = np.arange(np.sum(is_from_origin))
            output['Cluster_ID'] = clusters[is_from_origin]
            if args.write_pc:
                for i, pci in enumerate(pc):
                    output['PC%02i'%(i+1)] = pci[is_from_origin]
            output.to_csv(cl_file_name)

        # Write summary information
        summary_name = args.output_base 
        summary_name += '_n%02i_s%02i_k%02i_summary.csv'%(args.n_components, args.random_state, k)
        print('Writing summary information to', summary_name)
        output = pd.DataFrame()
        results = [cids, sizes, cc_orig_fn, cc_orig_sim, cc_orig_id]
        columns = ['Cluster_ID', 'Size', 'Centroid_Original_File', 'COrig_File_ID', 'Centroid_Original_Index']
        for col, res in zip(columns, results):
            output[col] = res
        output.to_csv(summary_name, index=False)
        print(output)

    # Write information about all k values in this study
    file_name_ssd = args.output_base 
    file_name_ssd += '_n%02i_s%02i_ssd.csv'%(args.n_components, args.random_state)
    print('Writing the sums of the squared distances to', file_name_ssd)
    ssd = pd.DataFrame()
    ssd['Num_Clusters'] = args.n_clusters
    ssd['Sum_Squ_Dist'] = sum_sqrd 
    ssd.to_csv(file_name_ssd, index=False)

    # TODO: Create plots or save

