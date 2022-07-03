#!/usr/bin/env python3

import pandas as pd
from schrodinger.application.desmond.packages import analysis
from trajectory_io import load_trajectory, get_an_aid


def get_an_aid(cms_model, asl_string):
    aid_list = cms_model.select_atom(asl_string)
    if (len(aid_list) != 1):
        raise Exception("The number of chosen atoms is %d."%(len(aid_list)))
    else:
        return(aid_list[0])


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
    msys_model, cms_model, tr = load_trajectory(args.cms_file, args.xtc_file, args.trj_dir)
                        
    # time
    frame_time = []
    for item in tr:
        frame_time.append(item.time)

    # define residue names
    if args.residue_names is None:
        args.residue_names = args.residue_numbers
    else:
        assert len(args.residue_names) == len(args.residue_numbers)

    # define analyzers
    analyzers = []
    distance_names = []
    for i in range(len(args.residue_numbers)):
        for j in range(i+1, len(args.residue_numbers)):
            rnum_i, rname_i  = args.residue_numbers[i], args.residue_names[i]
            rnum_j, rname_j  = args.residue_numbers[j], args.residue_names[j]
            first_asl = '(res.num %i) AND (atom.ptype " CA ") AND (chain.name %s)'%(rnum_i, args.chain_id)
            first_aid = get_an_aid(cms_model, first_asl)
            second_asl = '(res.num %i) AND (atom.ptype " CA ") AND (chain.name %s)'%(rnum_j, args.chain_id)
            second_aid = get_an_aid(cms_model, second_asl)
            analyzers.append(analysis.Distance(msys_model, cms_model, first_aid, second_aid))
            distance_names.append("%s-%s"%(rname_i, rname_j))
 
    #compute result
    dists = analysis.analyze(tr, *analyzers)

    # create and write csv file
    results = pd.DataFrame()
    results['time (ps)'] = frame_time
    if (len(distance_names) > 1):  # dists is a two dimensional array
        for name, dist_values in zip(distance_names, dists):
            results[name] = dist_values
    else:
        results[distance_names[0]] = dists
    results.to_csv(args.outputfile, index=False)
