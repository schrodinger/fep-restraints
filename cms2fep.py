"""
This script converts relaxed membrane cms system to a system ready for FEP
run.
"""
import argparse
import os

import schrodinger.structure as structure
from schrodinger.application.desmond import constants
from schrodinger.application.desmond.packages import cui
from schrodinger.structutils import analyze

KEEP_CT_TYPE = [
    constants.CT_TYPE.VAL.MEMBRANE, constants.CT_TYPE.VAL.SOLUTE,
    constants.CT_TYPE.VAL.SOLVENT
]

ION_CT_TYPE = [
    constants.CT_TYPE.VAL.ION, constants.CT_TYPE.VAL.POSITIVE_SALT,
    constants.CT_TYPE.VAL.NEGATIVE_SALT
]


def _extract_ligand(st, lig_asl):
    """
    Extract ligand if ASL is provided, and return receptor_ct and ligand_ct
    Delete ffio related properties.
    """
    if constants.CT_TYPE in st.property:
        del st.property[constants.CT_TYPE]

    solute_ct = st.extract(list(range(1, st.atom_total + 1)), copy_props=True)
    ligand_ct = None

    if lig_asl:
        lig_idx = analyze.evaluate_asl(st, lig_asl)
        if lig_idx:
            ligand_ct = st.extract(lig_idx, copy_props=True)
            solute_ct = st.extract(analyze.evaluate_asl(st, 'not ' + lig_asl),
                                   copy_props=True)
            print("\tFound ligand %s with %i atoms" %
                  (ligand_ct.atom[1].pdbres.strip(), ligand_ct.atom_total))
        else:
            print("\tLigand ASL failed to match any ligands.")
    return solute_ct, ligand_ct


def _cms_to_fep_inputs(st_list, lig_asl=None):
    """
    Extracts solute, water and membrane CTs.

    The return list contains cts in the following order:

    - receptor
    - membrane
    - solvent
    - ligand (if ligand_asl provided)

    """
    solvent_ct = None
    membrane_ct = None
    solute_ct = None
    ligand_ct = None
    for st in st_list:
        prop = st.property.get(constants.CT_TYPE)
        if prop == constants.CT_TYPE.VAL.MEMBRANE:
            membrane_ct = st
        elif prop == constants.CT_TYPE.VAL.SOLVENT:
            solvent_ct = st
        elif prop == constants.CT_TYPE.VAL.SOLUTE:
            solute_ct, ligand_ct = _extract_ligand(st, lig_asl)
    for st in st_list:
        if st.property.get(constants.CT_TYPE) in ION_CT_TYPE:
            solute_ct = solute_ct.merge(st, copy_props=False)
    solute_ct.title = 'receptor'
    if ligand_ct:
        ligand_ct.title = 'ligand-0'
    if solvent_ct:
        for atom in solvent_ct.atom:
            atom.property.pop(constants.FEP_ABSOLUTE_ENERGY, None)
    if membrane_ct is None:
        out_sys = [st for st in [solute_ct, solvent_ct, ligand_ct] if st]
    else:
        out_sys = [st for st in [solute_ct, membrane_ct, solvent_ct, ligand_ct] if st]
    return out_sys


def main(argv=None):
    parser = argparse.ArgumentParser(description="convert cms to FEP pv file.")
    parser.add_argument('cms_file',
                        metavar='<cms>',
                        help="CMS file of relaxed membrane system")
    parser.add_argument('-ligand',
                        dest='ligand_asl',
                        type=cui.validate_and_enclose_asl,
                        help="If the system contains a ligand that, and you "
                        "want to run FEP on that ligand, please specify "
                        "its ASL")
    parser.add_argument('-o',
                        dest='output_fn',
                        metavar='*_pv.mae',
                        action='store',
                        default='fep_pv.mae',
                        help='Output file to be used as FEP input.')
    args = parser.parse_args(argv)
    if not os.path.isfile(args.cms_file):
        cui.error("CMS file (%s) is not found." % args.cms_file)
    if not args.cms_file.endswith('.cms'):
        cui.error("Input file should be in CMS format.")

    if args.output_fn.endswith('.mae'):
        if not args.output_fn.endswith('_pv.mae'):
            args.output_fn = args.output_fn.replace('.mae', '_pv.mae')
    else:
        args.output_fn += '_pv.mae'

    streader = structure.StructureReader(args.cms_file)
    cms_components = [st for st in streader]
#    if constants.CT_TYPE.VAL.MEMBRANE not in [
#            st.property.get(constants.CT_TYPE) for st in cms_components
#    ]:
#        cui.error("Provided system does not contain membrane-CT.")

    output_cts = _cms_to_fep_inputs(cms_components, args.ligand_asl)
    stwriter = structure.StructureWriter(args.output_fn)
    for ct in output_cts:
        # remove link to the trajectory
        if 's_chorus_trajectory_file' in ct.property:
            del ct.property['s_chorus_trajectory_file']
        stwriter.append(ct)
    stwriter.close()
    print("Output %s has been written." % args.output_fn)


if __name__ == '__main__':
    main()
