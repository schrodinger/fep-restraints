#!/usr/bin/env python3

import os, os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from schrodinger.application.desmond.packages import traj, topo
from tasks.io_trajectory import write_frames, copy_topology, extract_frames_by_value
from tasks.io_features import write_features_to_csv, read_features_from_csv_files, sort_features, \
    calculate_ca_distances, calculate_backbone_torsions, calculate_sidechain_torsions
from tasks.comparison import plot_most_different_distributions, relative_entropy_analysis
from tasks.clustering_on_pca import kmeans_on_pca, scatterplot_pca_by_system, plot_pca_by_system, elbow_plot, pc_cluster_plot
from tasks.rmsf_from_trajectories import calculate_rmsf, write_coordinates, plot_cluster_rmsf, select_subset_model


def calculate_features(simulations, args, feature_type='ca-distance', chain_id_in_name=False):

    feature_files = []

    for i in simulations.index:

        # Get names and info for this simulation
        top_file = simulations['Topology'][i]
        trj_file = simulations['Trajectory'][i]
        sys_name = simulations['System_Name'][i]
        sys_prop = selections[selections['System_Name']==sys_name]
        chain_id = list(sys_prop['BindingPocket_Chain'])
        res_nums = [np.array(rn.split(' '),dtype=int) for rn in sys_prop['BindingPocket_ResNum']] 
        # Print names and info for this simulation
        print('Top. File:', top_file)
        print('Trj. File:', trj_file)
        print('System: '+sys_name)
        for _i, _chain in enumerate(chain_id):
            print('Chain:', _chain)
            print('Residues:', res_nums[_i])

        # Load the simulation data
        msys_model, cms_model = topo.read_cms(top_file)
        trj = traj.read_traj(trj_file) 
        
        # Calculate the distances between the C-alpha atoms
        if feature_type == 'ca-distance':
            time, feat_names, feat_values = calculate_ca_distances(
                msys_model, cms_model, trj, chain_id, res_nums,
                chain_id_in_name=chain_id_in_name
            )
        elif feature_type == 'bb-torsion':
            time, feat_names, feat_values = calculate_backbone_torsions(
                msys_model, cms_model, trj, chain_id, res_nums,
                chain_id_in_name=chain_id_in_name
            )
        elif feature_type == 'sc-torsion':
            time, feat_names, feat_values = calculate_sidechain_torsions(
                msys_model, cms_model, trj, chain_id, res_nums,
                chain_id_in_name=chain_id_in_name
            )       
        # ... and write them to a CSV file
        out_file = os.path.join(args.output_dir,f'1-features/%ss_%04i.csv' % (feature_type, i))
        feature_files.append(out_file)
        write_features_to_csv(out_file, feat_values, feat_names, time)
        print('Wrote %ss from %s to %s.\n'%(feature_type, trj_file, out_file) )

    return feature_files


def compare_features(simulations, args, feature_files_key, feature_type='ca-distance',
                     output_name='features', out_column='FeatureName',
                     sim_label_a='active', sim_label_b='inactive'):
    # Read the input files
    names, data, origin, orig_id = read_features_from_csv_files(simulations[feature_files_key])
    # Split data in active and inactive (according to their origin simulation)
    origin_a = list(simulations[simulations['Start_Label']==sim_label_a].index)
    origin_b = list(simulations[simulations['Start_Label']==sim_label_b].index)
    print('Simulations %s: ' % sim_label_a, origin_a)
    print('Simulations %s: ' % sim_label_b, origin_b)
    print('Shape of the total data:   ', data.shape)
    is_a = [o in origin_a for o in origin]
    is_b = [o in origin_b for o in origin]
    data_a = data[is_a]
    data_b = data[is_b]
    print('Shape of the %s data:' % sim_label_a, data_a.shape)
    print('Shape of the %s data:' % sim_label_b, data_b.shape) 
    # Run the relative-entropy analysis
    data_names, jsd, kld_ab, kld_ba = relative_entropy_analysis(
        names, names, data_a, data_b, 
        bin_width=None, bin_num=10, verbose=True
    )
    # Sort the features by how much their distributions differ
    jsd_sorted = sort_features(data_names, jsd)
    out_data = pd.DataFrame(jsd_sorted, columns=[out_column,'JSD'])
    out_csv = os.path.join(args.output_dir, f'2-comparison/%s_sorted-by-jsd.csv' % output_name)
    out_data.to_csv(out_csv)
    # Plot the 20 most different distributions
    out_plot = os.path.join(args.output_dir, f'2-comparison/%s_largest-jsd' % output_name)
    plot_most_different_distributions(
        jsd_sorted, names, names, data_a, data_b, out_plot,
        showstart=args.showstart, feature_type = feature_type
    )
    return jsd_sorted


if __name__ == "__main__":
    """ 
    Select cluster centroids, averages, and RMSF via PCA and k-means clustering on C-alpha distances between all selected residues. 
    """
       
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_file',  type=str, help='Name of the input file in csv format determining topology+trajectory files.')
    parser.add_argument('-s', dest='select_file', type=str, help='Name of the input file in csv format determining residue selections for binding pocket and restraints.') 
    parser.add_argument('-o', dest='output_dir',  type=str, help='Name of the output directory', default='clustering_results')       
    parser.add_argument('-n', dest='n_components',type=int, help='Number of top principal components to consider', default=3)     
    parser.add_argument('-k', nargs='+', dest='n_clusters', type=int, help='numbers k of clusters to attempt (arbitrary number)', default=[1,2,3,4]) 
    parser.add_argument('-r', dest='random_state', type=int, help='random state for k-means algorithm', default=42)    
    parser.add_argument('-w', dest='write_traj', action='store_true', default=False, help='Sort all trajectory frames by cluster. Warning: Currently this only works if all trajectories have the same topology and otherwise causes the code to fail.')
    parser.add_argument('-t', dest='threshold', type=float, default=None)
    parser.add_argument('--showstart', dest='showstart', action='store_true', default=False)
    parser.add_argument('--skip-comparison', dest='skip_comparison', action='store_true', default=False)
    parser.add_argument('--sim-label-a', dest='sim_label_a', type=str, default='active')
    parser.add_argument('--sim-label-b', dest='sim_label_b', type=str, default='inactive')
    parser.add_argument('--chain-id-in-name', dest='chain_id_in_name', action='store_true', default=False, help='Store the chain ID in the feature name. Note: For this to work, the chain naming has to be consistent across all input simulations!')
    args = parser.parse_args()

    # Read the input files
    print("\n* - Reading the input files. - *\n")
    simulations = pd.read_csv(args.input_file, comment='#')
    selections = pd.read_csv(args.select_file, comment='#')
    print(simulations)
    print(selections)

    # Create the output directories
    os.makedirs(args.output_dir, exist_ok=True)
    steps = ['1-features', '2-comparison', '3-pca', '4-clustering', '5-rmsf']
    if args.write_traj:
        steps += ['6-sorted']
    for subdir in steps:
        newdir = os.path.join(args.output_dir, subdir)
        os.makedirs(newdir, exist_ok=True)  

    # * ------------------------------------------------------------------ * #
    # *  Calculate the C-alpha distances and backbone/sidechain torsions.  * #
    # * ------------------------------------------------------------------ * #
  
    print("\n* - Calculating the features. - *\n")

    simulations['CA-Dist_File'] = calculate_features(
        simulations, args, 'ca-distance',
        chain_id_in_name=args.chain_id_in_name
    )
    simulations['BB-Tors_File'] = calculate_features(
        simulations, args, 'bb-torsion',
        chain_id_in_name=args.chain_id_in_name
    )
    simulations['SC-Tors_File'] = calculate_features(
        simulations, args, 'sc-torsion',
        chain_id_in_name=args.chain_id_in_name
    )
 
    # * ------------------------------------------ * #
    # *  Compare the features of the simulations.  * #
    # * ------------------------------------------ * #

    if not args.skip_comparison:

        print("\n* - Comparing the C-alpha distances of the simulations. - *\n")
        _ = compare_features(
            simulations, args, 'CA-Dist_File', feature_type='ca-distance',
            output_name='ca-dist', out_column='CA-Distance',
            sim_label_a=args.sim_label_a, sim_label_b=args.sim_label_b
        )
        
        print("\n* - Comparing the backbone torsions of the simulations. - *\n")
        _ = compare_features(
            simulations, args, 'BB-Tors_File', feature_type='bb-torsion',
            output_name='bb-tors', out_column='BB-Torsion',
            sim_label_a=args.sim_label_a, sim_label_b=args.sim_label_b
        )    
        
        print("\n* - Comparing the sidechain torsions of the simulations. - *\n")
        _ = compare_features(
            simulations, args, 'SC-Tors_File', feature_type='sc-torsion',
            output_name='sc-tors', out_column='SC-Torsion',
            sim_label_a=args.sim_label_a, sim_label_b=args.sim_label_b
        )

    # * ------------------------------ * #
    # *  Principal Component Analysis  * #
    # * ------------------------------ * #

    paramstr = 's%02i'%(args.random_state)

    print("\n* - Determining the principal components. - *\n")
    names, data, origin, orig_id = read_features_from_csv_files(simulations['CA-Dist_File'])
    pca = PCA(random_state=args.random_state)
    pca.fit(data.T)
    pc = pca.components_
    ev_ratio = pca.explained_variance_ratio_
    
    # Write the explained variance ratio to a CSV file and to stdout
    ev_csv = os.path.join(args.output_dir,'3-pca/ca-distances_pca_'+paramstr+'_ev-ratio.csv')
    ev_output = pd.DataFrame()
    ev_output['PC'] = 1 + np.arange(len(ev_ratio), dtype=int)
    ev_output['Explained_Variance_Ratio'] = ev_ratio
    ev_output.to_csv(ev_csv, index=False)
    print('Explained variance ratio:')
    for i, evr in enumerate(ev_ratio[:args.n_components]):
        print(' PC%02i: %1.4f'%(i+1, evr))

    # Write PCA results to CSV file
    out_csv = os.path.join(args.output_dir,'3-pca/ca-distances_pca_'+paramstr+'.csv')
    pca_output = pd.DataFrame()
    pca_output['Origin'] = origin
    max_pc = max([12, args.n_components])
    for i, pci in enumerate(pc[:max_pc]):
        pca_output['PC%i'%(i+1)] = pci
    pca_output.to_csv(out_csv, index=False) 

    # Plot PCA results by origin system
    out_pdf = os.path.join(args.output_dir,'3-pca/ca-distances_pca_'+paramstr)
    plot_pca_by_system(
        pc[:max_pc], origin, simulations, out_pdf, showstart=args.showstart
    )
    scatterplot_pca_by_system(
        pc, 0, 1, origin, simulations, out_pdf+'_pc1and2', showstart=args.showstart
    )
    scatterplot_pca_by_system(
        pc, 0, 2, origin, simulations, out_pdf+'_pc1and3', showstart=args.showstart
    )
    scatterplot_pca_by_system(
        pc, 1, 2, origin, simulations, out_pdf+'_pc2and3', showstart=args.showstart
    )

    # For all further work, only keep the requested number of PCs
    pc = pc[:args.n_components]
    # Update the parameter string to reflect this
    paramstr = 'n%02i_%s'%(args.n_components, paramstr)


    # * -------------------- * #
    # *  K-Means Clustering  * #
    # * -------------------- * #
    
    # Perform k-means clustering in PC space for various k values.
    sum_sqrd = []
    cl_files = []
    sum_file = []
    centroid_files = []
    cluster_centers = []
    for ik, k in enumerate(args.n_clusters):

        print('\n* - Running k-means clustering with k=%i - *\n'%k)

        paramstr_k = '%s_k%02i'%(paramstr, k)
        outputf = os.path.join(args.output_dir,'4-clustering/pca-kmeans_'+paramstr_k)
        
        # Run the k-means clustering
        cids, sizes, centers, cdist, cc_orig_sim, cc_orig_id, inertia, cl_files_k, sum_file_k = kmeans_on_pca(
            pc, k, args.random_state, origin, orig_id, output_base=outputf, input_files=None, write_pc=True
            )
        sum_sqrd.append(inertia)
        cl_files.append(cl_files_k)
        sum_file.append(sum_file_k)

        # Get the trajectory and save the frame with the centroid
        centroid_files_k = []
        for cl_id, c_file, c_frame in zip( cids, cc_orig_sim, cc_orig_id ):
            top_file = simulations['Topology'][c_file]
            trj_file = simulations['Trajectory'][c_file]
            _, top = topo.read_cms(top_file)
            trj = traj.read_traj(trj_file)
            out_dir = os.path.join(args.output_dir,'4-clustering')
            out_fn = 'pca-kmeans_'+paramstr_k+'_centroid%02i'%cl_id
            write_frames(top, trj, [c_frame], out_dir, frame_names=[out_fn])
            cf = os.path.join(out_dir, out_fn+'.cms')
            centroid_files_k.append(cf)
            print('Wrote centroid #%i to file %s'%(cl_id, cf))
        centroid_files.append(centroid_files_k)
        cluster_centers.append(centers)

    # Write information about all k values in this study
    file_name_ssd = os.path.join(args.output_dir,'4-clustering/pca-kmeans_'+paramstr+'_ssd.csv')
    print('Writing the sums of the squared distances to', file_name_ssd)
    ssd = pd.DataFrame()
    ssd['Num_Clusters'] = args.n_clusters
    ssd['Sum_Squ_Dist'] = sum_sqrd 
    ssd.to_csv(file_name_ssd, index=False)

    # Elbow plot for SSD over k
    plot_name_ssd = os.path.join(args.output_dir,'4-clustering/pca-kmeans_'+paramstr+'_ssd')
    elbow_plot(args.n_clusters, sum_sqrd, plot_name_ssd)

    # Plot clusters in PCA space
    for ik, k in enumerate(args.n_clusters):
        centers = cluster_centers[ik]
        cl_files_k = cl_files[ik]
        paramstr_k = '%s_k%02i'%(paramstr, k)
        out_name_k = '4-clustering/pca-kmeans_'+paramstr_k+'_pc-clusters'
        out_pca_cl = os.path.join(args.output_dir, out_name_k)
        pc_cluster_plot(pc, cl_files_k, centers, out_pca_cl)


    # * ------ * #
    # *  RMSF  * #
    # * ------ * # 

    # Go through all clustering runs
    rmsf_files = []
    for ik, k in enumerate(args.n_clusters):  

        rmsf_files_k = []
        paramstr_k = '%s_k%02i'%(paramstr, k)
        csv_files = cl_files[ik]
        top_files = simulations['Topology']
        trj_files = simulations['Trajectory']
        sys_names = simulations['System_Name']
        summary = pd.read_csv(sum_file[ik])
        centroid_origin_file_id = summary['COrig_File_ID']
        centroid_origin_sysname = sys_names[centroid_origin_file_id]
        centroid_files_k = centroid_files[ik]

        # Go through all cluster IDs 
        for value in range(k):

            print('\n* - Calculating RMSF for cluster %i from k-means with k=%i - *\n'%(value, k))            
            paramstr_cl = '%s_cluster%02i'%(paramstr_k, value)

            # Centroid and its properties
            centroid_file = centroid_files_k[value]
            cc_topol_file = top_files[centroid_origin_file_id[value]]
            cc_origin_sys = sys_names[centroid_origin_file_id[value]]
            # Construct selection strings
            cc_selections = selections[selections['System_Name']==cc_origin_sys]
            # ... for alignment 
            ref_asl_align = ''
            for _chain, _resnumstr in zip(cc_selections['BindingPocket_Chain'], 
                                          cc_selections['BindingPocket_ResNum']):
                ref_asl_align += '( chain name '+_chain+' AND res.num '+_resnumstr+' ) OR '
            ref_asl_align = ref_asl_align[:-4] # cut the final OR
            # ... and for writing restraints
            ref_asl_write = ''
            for _resasl in cc_selections['Restraints_Res_ASL']:
                ref_asl_write += '( '+_resasl+' ) OR '
            ref_asl_write = ref_asl_write[:-4] # cut the final OR 
            print(' - Reference -')
            print(' File: %s\n Align ASL: %s\n Write ASL: %s\n'%(centroid_file, ref_asl_align, ref_asl_write))

            print(' - Trajectories -')
            cluster_top_files, cluster_trj_files, cluster_csv_files = [], [], [] 
            cluster_asl_align, cluster_asl_write = [], []
            for top, trj, csv, sys in zip(top_files, trj_files, csv_files, sys_names):
                # Use only the simulations of the same system as the centroid
                if sys != cc_origin_sys:
                    continue
                cluster_top_files.append(top)
                cluster_trj_files.append(trj)
                cluster_csv_files.append(csv)
                # Construct selection strings
                cluster_selections = selections[selections['System_Name']==sys] 
                # ... for alignment
                asl_align = ''
                for _chain, _resnumstr in zip(cluster_selections['BindingPocket_Chain'], 
                                              cluster_selections['BindingPocket_ResNum']):
                    asl_align += '( chain name '+_chain+' AND res.num '+_resnumstr+' ) OR '
                asl_align = asl_align[:-4] # cut the final OR
                # ... and for writing restraints.
                asl_write = ''
                for _resasl in cluster_selections['Restraints_Res_ASL']:
                    asl_write += '( '+_resasl+' ) OR '
                asl_write = asl_write[:-4] # cut the final OR 
                # Append the selection strings.
                cluster_asl_align.append(asl_align)
                cluster_asl_write.append(asl_write)
                print(' Top. File: %s\n Trj. File: %s\n Align ASL: %s\n Write ASL: %s\n'%(top, trj, asl_align, asl_write))

            # Calculate the RMSF of this cluster
            rmsf_per_atom, pos_average, cms_model_ref_new = calculate_rmsf(
                centroid_file, cluster_top_files, cluster_trj_files, cluster_csv_files, value, 
                ref_asl_align, ref_asl_write, cluster_asl_align, cluster_asl_write, 
                align_avg=True, threshold=args.threshold)   

            # Write the RMSF to a CSV file.  
            output = pd.DataFrame()
            output['RMSF'] = rmsf_per_atom
            if args.chain_id_in_name:
                output['chain'] = [a.chain for a in cms_model_ref_new.atom]
            output['pdbres'] = [a.pdbres for a in cms_model_ref_new.atom]
            output['resnum'] = [a.resnum for a in cms_model_ref_new.atom]
            output['pdbname'] = [a.pdbname for a in cms_model_ref_new.atom]
            out_csv_file = os.path.join(args.output_dir,'5-rmsf/pca-kmeans_'+paramstr_cl+'_rmsf.csv')
            output.to_csv(out_csv_file, index=False)
            rmsf_files_k.append(out_csv_file)

            # Write the RMSF on the input topology structure.
            _, cms_model_top = topo.read_cms(cluster_top_files[0])
            aidlist_write_top = cms_model_top.select_atom(str(ref_asl_write))
            cms_model_top_new = select_subset_model(cms_model_top, aidlist_write_top)
            out_fn_top = os.path.join(args.output_dir,'5-rmsf/pca-kmeans_'+paramstr_cl+'_rmsf_top.cms')
            _ = write_coordinates(out_fn_top, cms_model_top_new, xyz=None, sigma=rmsf_per_atom)

            # Write the RMSF on the reference structure (the centroid).
            out_fn_ref = os.path.join(args.output_dir,'5-rmsf/pca-kmeans_'+paramstr_cl+'_rmsf_ref.cms')
            _ = write_coordinates(out_fn_ref, cms_model_ref_new, xyz=None, sigma=rmsf_per_atom)

            # Write the RMSF on the average structure.
            out_fn_avg = os.path.join(args.output_dir,'5-rmsf/pca-kmeans_'+paramstr_cl+'_rmsf_avg.cms')
            _ = write_coordinates(out_fn_avg, cms_model_ref_new, pos_average, sigma=rmsf_per_atom)

            # Write this cluster as a trajectory (xtc format)
            print('\n* - Writing cluster %i from k-means with k=%i as a trajectory - *\n'%(value, k))  
            if args.write_traj:
                cluster_output = os.path.join(args.output_dir,'6-sorted/pca-kmeans_'+paramstr_cl)  
                copy_topology(cc_topol_file, cluster_output+'.cms')
                extract_frames_by_value(cluster_trj_files, cluster_output+'.xtc', cluster_csv_files, value)
            print(' ')

        # Append the list of RMSF files
        rmsf_files.append(rmsf_files_k)

        # Plot every cluster of this clustering run
        # Only works if chain ID is neglected
        if not args.chain_id_in_name:
            out_name_k = '5-rmsf/pca-kmeans_'+paramstr_k+'_rmsf'
            out_rmsf_plot = os.path.join(args.output_dir, out_name_k)
            plot_cluster_rmsf(k, rmsf_files_k, names, out_rmsf_plot)

