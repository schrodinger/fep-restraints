from schrodinger.utils import sea
from schrodinger.application.desmond import cmj

from collections import defaultdict

def set_default(sea_map, key, default):
    '''Analagous to dict.setdefault.'''
    sea_map[key] = sea_map[key] if key in sea_map else default
    return sea_map[key]

def structure_to_restraints(asl: str, res_type: str, force_constant:float, sigma:float) -> dict:
    restraints = defaultdict(list)
    atom_restraint = {'atoms': asl,
                      'force_constants': force_constant
                     }    
    if (res_type == 'posre_harm'):
        restraints['posre_harm'].append(atom_restraint)        
    elif (res_type == 'posre_fbhw'):
        atom_restraint['sigma'] = sigma            
        restraints['posre_fbhw'].append(atom_restraint)
    else:
        pass
    return(restraints)

def add_restraints_to_stgs(stgs: sea.Map, restraints: dict):
    '''Add the restraints to the first assign_forcefield stage.

    We use the assign_forcefield stage so the terms will not be ignored
    or altered by `restraints.existing` parameters in later stages.
    '''
    for stg in stgs['stage']:
        # Only add once
        if stg.__NAME__ == 'assign_forcefield':
            restraint_map = set_default(stg, 'restraints', sea.Map())
            new_restraints = set_default(restraint_map, 'new', sea.List())
            for key, key_restraints in restraints.items():
                for restraint in key_restraints:
                    sea_restraint = sea.Map()
                    sea_restraint['name'] = sea.Atom(key)
                    for restraint_key, restraint_val in restraint.items():
                        sea_restraint[restraint_key] = sea.Atom(restraint_val)
                    new_restraints.append(sea_restraint)
            break

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', dest='msj_file', type=str, help='input msj file')
    parser.add_argument('-a', dest='atoms_asl', type=str, help='asl for the group of atoms', default='protein')
    parser.add_argument('-t', dest='res_type', type=str, help='restraint type: posre_harm; posre_fbhw')
    parser.add_argument('-k', dest='force_constant', type=float, help='force constant')
    parser.add_argument('-sigma', dest='sigma', type=float, help='sigma', default=None)    
    parser.add_argument('-o', dest='output_file', type=str, help='output structure file', default='result.msj')    
    args = parser.parse_args()    

    restraints = structure_to_restraints(args.atoms_asl, args.res_type, args.force_constant, args.sigma)

    stgs = cmj.msj2sea(args.msj_file) # read the msj file
    add_restraints_to_stgs(stgs, restraints)

    # Save the MSJ
    stgs.add_tag("setbyuser")
    cmj.write_sea2msj(stgs.stage, fname=args.output_file)
