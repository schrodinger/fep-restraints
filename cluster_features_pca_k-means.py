import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from features_io import read_features_from_csv_files

if __name__ == "__main__":
    """ 
    Perform k-means clustering in principal component space.
    """
        
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_files', nargs='+', type=str, help='input csv files')
    parser.add_argument('-o', dest='output_file', type=str, help='name of the output csv file', default='clusters.csv')            
    args = parser.parse_args()

    names, data, origin = read_features_from_csv_files(args.input_files) 

    print(data.shape)