#!/usr/bin/env python3

import pandas as pd
from schrodinger.application.desmond.packages import analysis
from io_trajectory import load_trajectory, get_an_aid
from io_features import write_features_to_csv


def get_an_aid(cms_model, asl_string):
    aid_list = cms_model.select_atom(asl_string)
    if (len(aid_list) != 1):
        raise Exception("The number of chosen atoms is %d."%(len(aid_list)))
    else:
        return(aid_list[0])

def calculate_ca_distances(tr, chain_id, residue_numbers, residue_names=None):
    # time
    frame_time = []
    for item in tr:
        frame_time.append(item.time)
    # define residue names
    if residue_names is None:
        residue_names = residue_numbers
    else:
        assert len(residue_names) == len(residue_numbers)
    # define analyzers
    analyzers = []
    distance_names = []
    for i in range(len(residue_numbers)):
        for j in range(i+1, len(residue_numbers)):
            rnum_i, rname_i  = residue_numbers[i], residue_names[i]
            rnum_j, rname_j  = residue_numbers[j], residue_names[j]
            first_asl = '(res.num %i) AND (atom.ptype " CA ") AND (chain.name %s)'%(rnum_i, chain_id)
            first_aid = get_an_aid(cms_model, first_asl)
            second_asl = '(res.num %i) AND (atom.ptype " CA ") AND (chain.name %s)'%(rnum_j, chain_id)
            second_aid = get_an_aid(cms_model, second_asl)
            analyzers.append(analysis.Distance(msys_model, cms_model, first_aid, second_aid))
            distance_names.append("%s-%s"%(rname_i, rname_j))
    #compute result
    distances = analysis.analyze(tr, *analyzers)
    return frame_time, distance_names, distances


if __name__ == "__main__":
    """ 
    Calculates the C-alpha distancses between all selected residues. 
    """
       
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='cms_file', type=str, help='cms file')
    parser.add_argument('-x', dest='xtc_file', type=str, help='xtc file', default=None)
    parser.add_argument('-t', dest='trj_dir', type=str, help='trajecotry directory', default=None)    
    parser.add_argument('-n', nargs='+', dest='residue_names', type=str, help='residue names', default=None)
    parser.add_argument('-r', nargs='+', dest='residue_numbers', type=int, help='residue numbers', default=None)
    parser.add_argument('-l', dest='chain_id', type=str, help='chain identifier', default='A')
    parser.add_argument('-o', dest='outputfile', type=str, help='name of the output csv file', default='results.csv')            
    args = parser.parse_args()

    # load data
    msys_model, cms_model, trajectory = load_trajectory(
        args.cms_file, 
        args.xtc_file, 
        args.trj_dir
        )
    # calculate the distances between the C-alpha atoms
    time, dist_names, distances = calculate_ca_distances(
        trajectory, 
        args.chain_id, 
        args.residue_numbers, 
        residue_names=args.residue_names
        )
    # write the distances to a CSV file
    write_features_to_csv(
        args.outputfile, 
        distances, 
        dist_names, 
        time
        )
