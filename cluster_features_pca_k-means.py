import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from torch import _sobol_engine_initialize_state_
from features_io import read_features_from_csv_files

if __name__ == "__main__":
    """ 
    Perform k-means clustering in principal component space.
    """
        
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_files', nargs='+', type=str, help='input csv files')
    parser.add_argument('-o', dest='output_base', type=str, help='base name of the output files', default='pca-clusters')
    parser.add_argument('-n', dest='n_components', type=int, help='number of top principal components to consider', default=3) 
    parser.add_argument('-k', dest='n_clusters', type=int, help='numbers k of clusters to attempt (arbitrary number)', default=[1,2,3,4])     
    parser.add_argument('-s', dest='random_state', type=int, help='random state for k-means algorithm', default=42)      
    args = parser.parse_args()

    names, data, origin, orig_id = read_features_from_csv_files(args.input_files)

    # Print some information
    print("Origin:")
    for o in np.unique(origin):
        n_frames = np.sum(origin==o)
        top_file = args.input_files[o]
        print(" - %s: %5i frames"%(top_file, n_frames))

    # Perform PCA on the data.
    pca = PCA(n_components=args.n_components)
    pca.fit(data.T)
    pc = pca.components_
    print(pca.explained_variance_ratio_) 

    # Perform k-means clustering in PC space.
    # Try various k values.
    sum_sqrd = []
    for ik, k in enumerate(args.n_clusters):
        kmeans = KMeans(
            n_clusters = k, 
            init = 'k-means++', 
            random_state = args.random_state
        )
        # Fit and transform the data to cluster-distance space.
        pc_new = kmeans.fit_transform(pc.T)
        # Calculate the sum of squared distances for this k.
        sum_sqrd.append(kmeans.inertia_) 
        # Determine clusters and their centers.
        clusters = kmeans.fit_predict(pc.T)
        centers = kmeans.cluster_centers_
        sizes = [np.sum(clusters==ic) for ic in range(k)]
        # Find centroids and their index and origin.
        centroids = np.argmin(pc_new, axis=0)
        cc_orig_sim = origin[centroids]
        cc_orig_idx = orig_id[centroids]
        # Write information about these clusters to a csv file
        file_name_k = args.output_base + '_n%02i_s%02i_k%02i.csv'%(args.n_components, args.random_state, k)
        print('Writing cluster information to', file_name_k)
        results = pd.DataFrame()
        results['Cluster_ID'] = np.arange(k)
        results['Size'] = sizes
        results['Centroid_Original_File'] = [args.input_files[o] for o in cc_orig_sim]
        results['Centroid_Original_Index'] = cc_orig_idx
        results.to_csv(file_name_k, index=False)
        print(cc_orig_sim, cc_orig_idx)

    # Write information about all k values in this study
    file_name_ssd = args.output_base + '_n%02i_s%02i_ssd.csv'%(args.n_components, args.random_state)
    print('Writing the sums of the squared distances to', file_name_ssd)
    ssd = pd.DataFrame()
    ssd['Num_Clusters'] = args.n_clusters
    ssd['Sum_Squ_Dist'] = sum_sqrd 
    ssd.to_csv(file_name_ssd, index=False)

    # TODO: Create plots or save necessary detailed information

