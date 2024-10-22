import argparse
import os
import sys

from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils import analyze, rmsd
from schrodinger.application.desmond import cmj

from tasks.restraints import add_restraints_to_stgs, structure_to_restraints

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
    parser.add_argument("-s", "--sf",
            default=1.0,
            type=float,
            help="Scaling factor for the sigma from the restraints structure. The half-width of the flat bottom will be sf*sigma.")
    parser.add_argument("-b", "--bw",
            default=None,
            type=float,
            help="Half width of the flat bottom (A). Overrides use of sigma from the restraints structure.")
    parser.add_argument("-r", "--reference",
            default=None,
            type=str,
            help="Reference structure to align the restraints structure to.")
    parser.add_argument("-w", "--write_restraints_structure",
            default=False,
            action='store_true',
            help="Write the (aligned if reference provided) structure with the restraints to an MAE file.")
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

    # Load the structure with the restraints
    st = StructureReader.read(args.structure)
    at = analyze.evaluate_asl(st, args.asl)
    # Align the structure to the reference if given
    if args.reference is not None:
        st_fixed = StructureReader.read(args.reference)
        at_fixed = analyze.evaluate_asl(st_fixed, args.asl)
        rmsd.superimpose(st_fixed, at_fixed, st, at)
    # Write the (aligned) restraints structure
    if args.write_restraints_structure:
        base, ext = os.path.splitext(args.out)
        st_out = f"{base}_restraints.mae"
        with StructureWriter(st_out) as writer:
            writer.append(st)

    # Create restraints from the structure
    restraints = structure_to_restraints(st, args.asl, args.fc, args.bw, args.sf)
    stgs = cmj.msj2sea(args.msj)
    add_restraints_to_stgs(stgs, restraints)
    # Save the MSJ
    stgs.add_tag("setbyuser")
    cmj.write_sea2msj(stgs.stage, fname=args.out)

if __name__ == "__main__":
    main()
    