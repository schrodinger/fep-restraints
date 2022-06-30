#!/usr/bin/env python3

def get_an_aid(cms_model, asl_string):
    aid_list = cms_model.select_atom(asl_string)
    if (len(aid_list) != 1):
        raise Exception("The number of chosen atoms is %d."%(len(aid_list)))
    else:
        return(aid_list[0])

if __name__ == "__main__":
    """ to calculate the distancses between one or multiple pairs of atoms from a trjactory 
    input file example:
    ----------------------------------------------------------------------------------------------------------------------------------------
    distance name,first asl,second asl
    asp121 and ser211 (A), (res.num 121) AND (atom.ptype " CA ") AND (chain.name A), (res.num 211) AND (atom.ptype " CA ") AND (chain.name A)
    asp121 and asn310 (A), (res.num 121) AND (atom.ptype " CA ") AND (chain.name A), (res.num 310) AND (atom.ptype " CA ") AND (chain.name A)
    asn329 and ser211 (A), (res.num 329) AND (atom.ptype " CA ") AND (chain.name A), (res.num 211) AND (atom.ptype " CA ") AND (chain.name A)
    asn329 and asn310 (A), (res.num 329) AND (atom.ptype " CA ") AND (chain.name A), (res.num 310) AND (atom.ptype " CA ") AND (chain.name A)
    ----------------------------------------------------------------------------------------------------------------------------------------
    """
    
    from schrodinger.application.desmond.packages import analysis, traj, topo
    import pandas as pd
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='cms_file', type=str, help='cms file')
    parser.add_argument('-x', dest='xtc_file', type=str, help='xtc file', default=None)
    parser.add_argument('-t', dest='trj_dir', type=str, help='trajecotry directory', default=None)    
    parser.add_argument('-n', dest='distance_name', type=str, help='distance name', default='distance (A)')
    parser.add_argument('-f', dest='first_asl', type=str, help='asl for the first atom')
    parser.add_argument('-s', dest='second_asl', type=str, help='asl for the second atom')
    parser.add_argument('-i', dest='input_file', type=str, help='input csv file', default=None)
    parser.add_argument('-o', dest='outputfile', type=str, help='name of the output csv file', default='results.csv')            
    args = parser.parse_args()

    # load data
    msys_model, cms_model = topo.read_cms(args.cms_file)
    if ((not args.xtc_file) and (not args.trj_dir)):
        raise Exception("Neither the xtc file nor the trajectory directory is found.")
    elif (args.xtc_file and args.trj_dir):
        raise Exception("Both the xtc file and the trajectory directory are found.")
    elif (args.xtc_file):
        tr = traj.read_traj(args.xtc_file)  # xtc format trajectory
    else:
        tr = traj.read_traj(args.trj_dir)  # trj format trajectory
                        
    # time
    frame_time = []
    for item in tr:
        frame_time.append(item.time)

    # define analyzers
    analyzers = []
    distance_names = []
    if (args.input_file is None):
        first_aid = get_an_aid(cms_model, args.first_asl)
        second_aid = get_an_aid(cms_model, args.second_asl)
        analyzers.append(analysis.Distance(msys_model, cms_model, first_aid, second_aid))
        distance_names.append(args.distance_name)
    else:
        input_data = pd.read_csv(args.input_file)
        for pd_index, pd_row in input_data.iterrows():
            first_aid = get_an_aid(cms_model, pd_row['first asl'])
            second_aid = get_an_aid(cms_model, pd_row['second asl'])
            analyzers.append(analysis.Distance(msys_model, cms_model, first_aid, second_aid))
            distance_names.append(pd_row['distance name'])            
    
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
    results.to_csv(args.outputfile)
