#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os.path

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='def_file', type=str, help='csv file with centroid definition')
    parser.add_argument('-c', dest='cluster_id', type=int, help='cluster index')    
    args = parser.parse_args()

    def_df = pd.read_csv(args.def_file)
    origin = def_df['Centroid_Original_File']
    origin = np.array(origin)[args.cluster_id]
    system = os.path.basename(origin)
    print(system)