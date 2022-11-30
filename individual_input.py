import os
from pathlib import Path
import schrodinger.structure as structure
from schrodinger.application.desmond import constants
from schrodinger.application.desmond.starter.generator.abfep import prepare_inputs, generate

KEEP_CT_TYPE = [
    constants.CT_TYPE.VAL.MEMBRANE, 
    constants.CT_TYPE.VAL.SOLVENT
]

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Extract individual pv files from one FEP pv file.")
    parser.add_argument('mae_file', help="File with input structures for FEP")
    parser.add_argument('dir_base', help="Name for the output directory")
    args = parser.parse_args()

    # Read the file with input structures
    with structure.StructureReader(args.mae_file) as streader:
        cms_components = [st for st in streader]

    # Start the list of template structures. 
    # The first entry is supposed to be the receptor
    template = [cms_components[0]]
    # Start the list of ligands
    ligands = []
    # Add membrane and solvent to template.
    # Everything else is interpreted as a ligand.
    for st in cms_components[1:]:
        prop = st.property.get(constants.CT_TYPE)
        if prop in KEEP_CT_TYPE:
            template.append(st)
        else:
            ligands.append(st)

    # for each ligand
    for lig_st in ligands:
        # Get its name
        name = lig_st.title.replace(' ','_')
        # Create the job dir and go there
        lig_path = Path(args.dir_base, name)
        os.makedirs(lig_path, exist_ok=True)
        os.chdir(lig_path)
        # Write an input structure file
        mae_path = Path(name+'.mae')
        with structure.StructureWriter(mae_path) as stwriter:
            for st in template:
                stwriter.append(st)
            stwriter.append(lig_st)
        # Go back up in the dir hierarchy
        os.chdir('../..')
