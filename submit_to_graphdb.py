import argparse
import os
import sys
import numpy as np

from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils import analyze, rmsd
from schrodinger.application.desmond import cmj

from tasks.restraints import structure_to_restraints, submit_graphdb_job_with_restraints, check_distances


def parse_cmdline(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("pv_or_fmp_file",
            help="Filename for the PV structure or FMP file to submit to GraphDB.")
    parser.add_argument("yaml_file",
            help="Filename for the YAML file with the parameters for the GraphDB job.")
    parser.add_argument("restraint_file",
            help="Filename for the restraint structure. Coordinates taken as restraint centers and RMSF encoded in atom.r_desmond_sigma")
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
            help="Reference structure to align the restraints structure to. If none is provided, the PV or FMP file will be used.")
    parser.add_argument("-na", "--no_align",
            dest='align',
            default=True,
            action='store_false',
            help="Do not align the restraints structure to the reference.")
    parser.add_argument("-w", "--write_restraints_structure",
            default=False,
            action='store_true',
            help="Write the (aligned if reference provided) structure with the restraints to an MAE file.")
    parser.add_argument("--out",
            default=None,
            help="MSJ output name")
    args = parser.parse_args(argv)

    if args.out is None:
        base, ext = os.path.splitext(args.msj)
        args.out = f"{base}_with_restraints{ext}"
    if args.reference is None:
        args.reference = args.pv_or_fmp_file
    return args


def main():

    args = parse_cmdline(sys.argv[1:])

    # Load the structure with the restraints and the structure to align them to
    st = StructureReader.read(args.restraint_file)
    at = analyze.evaluate_asl(st, args.asl)
    st_fixed = StructureReader.read(args.reference)
    at_fixed = analyze.evaluate_asl(st_fixed, args.asl)

    # Align the restraints (if requested) and calculate RMSD
    if args.align:
        print(f"Aligning the restraints structure to the starting structure.")
        restraints_rmsd = rmsd.superimpose(st_fixed, at_fixed, st, at)
    else:
        print("Not aligning the restraints structure to the starting structure.")
        restraints_rmsd = rmsd.calculate_in_place_rmsd(st_fixed, at_fixed, st, at)
    print(f"The RMSD between the starting structure and the restraint centers is {restraints_rmsd:.2f} A.")

    # Check distances between restraint centers and starting structure
    check_distances(st, st_fixed, at, at_fixed, sf=args.sf)

    # Write the (aligned or non-aligned) restraints structure
    if args.write_restraints_structure:
        base, ext = os.path.splitext(args.out)
        st_out = f"{base}_restraints.mae"
        with StructureWriter(st_out) as writer:
            writer.append(st)

    # Create restraints from the structure
    restraints = structure_to_restraints(st, args.asl, args.fc, args.bw, args.sf)
    # Submit the job to GraphDB with the poses from a FMP/PV file and the usual parameters in a yaml file
    submit_graphdb_job_with_restraints(args.fmp_or_pv_file, args.yaml_file, restraints)

if __name__ == "__main__":
    main()
    