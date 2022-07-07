#!/usr/bin/env python3

if __name__ == "__main__":
    """ to calculate the RMSF for the atoms selected """

    from schrodinger.application.desmond.packages import analysis, traj, topo, structure
    import pandas as pd

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='cms_file', type=str, help='cms file')
    parser.add_argument('-x', dest='xtc_file', type=str, help='xtc file', default=None)
    parser.add_argument('-t', dest='trj_dir', type=str, help='trajecotry directory', default=None)
    parser.add_argument('-s', dest='asl', type=str, help='asl for the atoms to calculate the RMSF', default='atom.ptype CA')
    parser.add_argument('-f', dest='fit_asl', type=str, help='asl for the atoms used to fit', default='backbone and not (atom.ele H)')
    parser.add_argument('-r', dest='ref_pos', type=str, help='reference posistion from a trajectory frame or a mae file', default='0')
    parser.add_argument('-o', dest='outputfile', type=str, help='name of the output csv file', default='results.csv')                
    args = parser.parse_args()

    # load trajectory
    msys_model, cms_model = topo.read_cms(args.cms_file)
    if ((not args.xtc_file) and (not args.trj_dir)):
        raise Exception("Neither the xtc file nor the trajectory directory is found.")
    elif (args.xtc_file and args.trj_dir):
        raise Exception("Both the xtc file and the trajectory directory are found.")
    elif (args.xtc_file):
        tr = traj.read_traj(args.xtc_file)  # xtc format trajectory
    else:
        tr = traj.read_traj(args.trj_dir)  # trj format trajectory

    # get aid and gid
    aids = cms_model.select_atom(args.asl)
    fit_aids = cms_model.select_atom(args.fit_asl)
    fit_gids = topo.aids2gids(cms_model, fit_aids, include_pseudoatoms=False)
        
    # load pos of the reference structure
    try:
        iframe = int(args.ref_pos)
        ref_pos = tr[iframe].pos(fit_gids)
    except ValueError:
        pass

    # # get residue id for each atom
    pdb_atom_name = []
    residue_number = []
    residue_name = []
    
    for indexi in aids:
        pdb_atom_name.append(cms_model.atom[indexi].property['s_m_pdb_atom_name'])
        residue_number.append(cms_model.atom[indexi].property['i_m_residue_number'])
        residue_name.append(cms_model.atom[indexi].property['s_m_pdb_residue_name'])
        # print(f"{cms_model.atom[indexi].property['s_m_pdb_atom_name']} {cms_model.atom[indexi].property['i_m_residue_number']} {cms_model.atom[indexi].property['s_m_pdb_residue_name']}")

    results = pd.DataFrame(list(zip(pdb_atom_name, residue_number, residue_name)), columns = ['pdb atom name', 'residue number', 'residue name'])    
        
    # analysis
    one_analysis = analysis.RMSF(msys_model, cms_model, aids, fit_aids, ref_pos)
    results['rmsf'] = analysis.analyze(tr, one_analysis)

    # write results
    results.to_csv(args.outputfile)
        
