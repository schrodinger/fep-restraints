#!/usr/bin/env python3

from schrodinger.application.desmond.packages import analysis, traj, topo
from schrodinger import structure
import numpy as np


def write_coordinates(out_fname, model, xyz=None, sigma=None):
    print("Writing to model with %i atoms."%(len(model.atom))) 
    # Update the coordinates. 
    if xyz is not None:
        model.fsys_ct.setXYZ(xyz)
    # Add the RMSF as the atom sigma
    if sigma is not None:
        for a, atom in enumerate(model.atom):
            atom.property['r_desmond_sigma'] = sigma[a]
    # Write the structure to the designated name.
    with structure.StructureWriter(out_fname) as writer:
        writer.append(model.fsys_ct)


if __name__ == "__main__":
    """ 
    Read sigma information from one cms file and write it to the structure from another one.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_structure', type=str, help='cms file')
    parser.add_argument('-s', dest='sigma_structure', type=str, help='cms file with sigma information', default=None)
    parser.add_argument('-o', dest='output_filename', type=str, help='name of the output file', default='.')        
    args = parser.parse_args()

    # Read the file with the structure
    _, input_cms_model = topo.read_cms(args.input_structure)
    # Read the file with sigma information
    _, sigma_cms_model = topo.read_cms(args.sigma_structure)
    # Make sure both models have the same number of atoms
    print(len(sigma_cms_model.atom), len(input_cms_model.atom))
    assert len(sigma_cms_model.atom) == len(input_cms_model.atom)
    # Extract sigma information from the CMS model
    sigma = np.empty(len(sigma_cms_model.atom))
    for a, atom in enumerate(sigma_cms_model.atom):
        sigma[a] = atom.property['r_desmond_sigma']
    # Write the RMSF on reference structure.
    _ = write_coordinates(args.output_filename, input_cms_model, xyz=None, sigma=sigma)
