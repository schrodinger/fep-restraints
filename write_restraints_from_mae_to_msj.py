import argparse
import os
import sys

from collections import defaultdict

from schrodinger.structure import Structure, StructureReader, _StructureAtom
from schrodinger.structutils import analyze
from schrodinger.utils import sea
from schrodinger.application.desmond import cmj

def construct_atom_asl(st: Structure, atom: _StructureAtom, asl: str) -> str:
    '''Returns an ASL that (should) uniquely identify the provided atom.
    '''
    # First try a simplified ASL
    simple_asl = f"res.ptype {atom.pdbres.strip()} AND res.num {atom.resnum} AND atom.ptype {atom.pdbname.strip()}"
    aids = analyze.evaluate_asl(st, simple_asl)
    if len(aids) == 1:
        return simple_asl
    # Then tack on to the original ASL
    suffix_asl = f"({asl}) AND ({simple_asl})"
    aids = analyze.evaluate_asl(st, suffix_asl)
    if len(aids) == 1:
        return suffix_asl
    raise UserWarning("Unable to find an appropriate ASL to uniquely identify an atom")

def structure_to_restraints(st: Structure, asl: str, fc=50.0) -> dict:
    '''Take a structure and create position restraint (posre) terms from it.

    If the matched atom has a r_desmond_sigma property it will be assigned to posre_fbhw, otherwise posre_harm.
    The coordinates of the atom are taken as the center of the restraint.
    '''
    restraints = defaultdict(list)
    aids = analyze.evaluate_asl(st, asl)
    print(f"ASL '{asl}' matched {len(aids)} atoms.")
    for aid in aids:
        atom = st.atom[aid]
        atom_restraint = {'atoms': construct_atom_asl(st, atom, asl),
                          'ref': (atom.x, atom.y, atom.z),
                          'force_constants': fc,
                         }
        # The property 'dictionary' has no equivalent to the get method to make this easier
        if ('r_desmond_sigma' in atom.property) and (atom.property['r_desmond_sigma']):
            atom_restraint['sigma'] = atom.property['r_desmond_sigma']
            restraints['posre_fbhw'].append(atom_restraint)
        else:
            restraints['posre_harm'].append(atom_restraint)
    return restraints

def set_default(sea_map, key, default):
    '''Analagous to dict.setdefault.'''
    sea_map[key] = sea_map[key] if key in sea_map else default
    return sea_map[key]

def add_restraints_to_stgs(stgs: sea.Map, restraints: dict):
    '''Add the restraints to the first assign_forcefield stage.

    We use the assign_forcefield stage so the terms will not be ignored
    or altered by `restraints.existing` parameters in later stages.
    '''
    for stg in stgs['stage']:
        # Only add once
        if stg.__NAME__ == "assign_forcefield":
            restraint_map = set_default(stg, 'restraints', sea.Map())
            new_restraints = set_default(restraint_map, 'new', sea.List())
            for key, key_restraints in restraints.items():
                for restraint in key_restraints:
                    sea_restraint = sea.Map()
                    sea_restraint['name'] = sea.Atom(key)
                    for restraint_key, restraint_val in restraint.items():
                        if restraint_key == 'ref':
                            sea_val = sea.List()
                            for v in restraint_val:
                                sea_val.append(sea.Atom(v))
                        else:
                            sea_val = sea.Atom(restraint_val)
                        sea_restraint[restraint_key] = sea_val
                    new_restraints.append(sea_restraint)
            break

def parse_cmdline(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("structure",
            help="Filename for the structure.  Coordinates taken as restraint centers and RMSF encoded in atom.r_desmond_sigma")
    parser.add_argument("msj",
            help="MSJ file to modify")
    parser.add_argument("-a", "--asl",
            default="all",
            help="ASL to define which atoms are restrained (default all)")
    parser.add_argument("-f", "--fc",
            default=50.0,
            type=float,
            help="Force constant (kcal/mol/A**2)")
    parser.add_argument("out",
            default=None,
            help="MSJ output name")
    args = parser.parse_args(argv)

    if args.out is None:
        base, ext = os.path.splitext(args.msj)
        args.out = f"{base}_with_restraints{ext}"
    return args

def main():
    args = parse_cmdline(sys.argv[1:])

    st = StructureReader.read(args.structure)
    restraints = structure_to_restraints(st, args.asl, args.fc)

    stgs = cmj.msj2sea(args.msj)
    add_restraints_to_stgs(stgs, restraints)

    # Save the MSJ
    stgs.add_tag("setbyuser")
    cmj.write_sea2msj(stgs.stage, fname=args.out)

if __name__ == "__main__":
    main()
    