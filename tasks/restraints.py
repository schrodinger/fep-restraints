import numpy as np
from collections import defaultdict

from schrodinger.structure import Structure, _StructureAtom
from schrodinger.structutils import analyze
from schrodinger.utils import sea
from schrodinger.application.desmond import cmj
from schrodinger.application.scisol.graphdb_interface import jws, driver
from schrodinger.application.desmond.multisim import parse
from schrodinger.application.desmond.packages import topo

import base64


def construct_atom_asl(st: Structure, atom: _StructureAtom, asl: str) -> str:
    """
        Returns an ASL that (should) uniquely identify the provided atom.

        Parameters:
        ----------
        st : schrodinger.structure.Structure
            The structure object.
        atom : schrodinger.structure._StructureAtom
            The atom object.
        asl : str
            The original ASL to match the atoms.

        Returns:
        -------
        str
            The ASL that uniquely identifies the atom.
    """
    # First try a simplified ASL
    simple_asl = f"chain.name {atom.chain.strip()} AND res.ptype {atom.pdbres.strip()} AND res.num {atom.resnum} AND atom.ptype {atom.pdbname.strip()}"
    aids = analyze.evaluate_asl(st, simple_asl)
    if len(aids) == 1:
        return simple_asl
    # Then tack on to the original ASL
    suffix_asl = f"({asl}) AND ({simple_asl})"
    aids = analyze.evaluate_asl(st, suffix_asl)
    if len(aids) == 1:
        return suffix_asl
    raise UserWarning("Unable to find an appropriate ASL to uniquely identify an atom")

def structure_to_restraints(st: Structure, asl: str, fc=50.0, bw=None, sf=1.0) -> dict:
    """
        Take a structure and create position restraint (posre) terms from it.
        If the matched atom has a r_desmond_sigma property it will be assigned to posre_fbhw, otherwise posre_harm.
        The coordinates of the atom are taken as the center of the restraint.

        Parameters:
        ----------
        st : schrodinger.structure.Structure
            The structure object.
        asl : str
            The ASL to match the atoms.
        fc : float, optional
            The force constant of the restraints. Default is 50.0.
        bw : float, optional
            The bandwidth of the restraints. If None, the sigma value of the atoms will be used. Default is None.
        sf : float, optional
            A scaling factor for the sigma values. Default is 1.0.

        Returns:
        -------
        dict
            A dictionary containing the restraints to be added to the MSJs. The keys are the type of 
            restraints (e.g. 'posre_fbhw', 'posre_harm') and the values are lists of dictionaries 
            containing the restraint information. Each dictionary contains the following keys: 'atoms', 
            'ref', 'force_constants', and 'sigma' (optional
    """
    restraints = defaultdict(list)
    aids = analyze.evaluate_asl(st, asl)
    print(f"ASL '{asl}' matched {len(aids)} atoms.")
    for aid in aids:
        atom = st.atom[aid]
        atom_restraint = {
            'atoms': construct_atom_asl(st, atom, asl),
            'ref': (atom.x, atom.y, atom.z),
            'force_constants': fc,
            }
        if bw is None: # read sigma from the structure
            # The property 'dictionary' has no equivalent to the get method to make this easier
            if ('r_desmond_sigma' in atom.property) and (atom.property['r_desmond_sigma']):
                atom_restraint['sigma'] = sf*atom.property['r_desmond_sigma']
                restraints['posre_fbhw'].append(atom_restraint)
            else:
                restraints['posre_harm'].append(atom_restraint)
        else: # use constant sigma
            atom_restraint['sigma'] = bw
            restraints['posre_fbhw'].append(atom_restraint)
    return restraints

def read_posre_fbhw(sigma_structure_file: str) -> dict:
    """
        Read the sigma values from a Desmond structure file and create restraints.

        Parameters:
        ----------
        sigma_structure_file : str
            The path to the Desmond structure file containing the sigma values.

        Returns:
        -------
        dict
            A dictionary containing the restraints to be added to the MSJs. The keys are the type of 
            restraints (e.g. 'posre_fbhw', 'posre_harm') and the values are lists of dictionaries 
            containing the restraint information. Each dictionary contains the following keys: 'atoms', 
            'ref', 'force_constants', and 'sigma' (optional
    """
    _, sigma_cms_model = topo.read_cms(sigma_structure_file)
    restraints = {'posre_fbhw':[]}
    for a, atom in enumerate(sigma_cms_model.atom):
        sigma = atom.property['r_desmond_sigma']
        if atom.pdbname.replace(" ", "") == 'CA':
            atom_restraint = {
                'atoms': f'at.ptype CA and res. {atom.resnum}',
                'ref': (atom.x, atom.y, atom.z),
                'force_constants': 50,
                'sigma': sigma
                }
            restraints['posre_fbhw'].append(atom_restraint)
    return restraints

def check_distances(st, st_fixed, at, at_fixed, sf=1.0):
    """
        Check the distances between the restraints and the atoms in the reference structure.

        Parameters:
        ----------
        st : schrodinger.structure.Structure
            The structure object of the starting coordinates.
        st_fixed : schrodinger.structure.Structure
            The structure object of the reference coordinates.
        at : list
            A list of atom indices in the starting coordinates.
        at_fixed : list
            A list of atom indices in the reference coordinates.
        sf : float, optional
            A scaling factor for the sigma values. Default is 1.0.

        Returns:
        -------
        None
    """
    assert len(at) == len(at_fixed)
    counter = 0
    output_string = 'Residue:  dist., sigma, diff.\n'
    for i,j in zip(at, at_fixed):
        atom_i = st.atom[i]
        atom_j = st_fixed.atom[j]
        xyz_i = np.array(atom_i.xyz)
        xyz_j = np.array(atom_j.xyz)
        sigma = sf*atom_i.property['r_desmond_sigma']
        d = np.sqrt(np.sum(xyz_j - xyz_i)**2)
        if d > sigma:
            output_string += f'{atom_i.pdbres}{atom_i.resnum}: {d:.3f}, {sigma:.3f}, {d-sigma:.3f}\n'
            counter += 1
    if counter == 0:
        print("All restraint centers are within the sigma range of the starting coordinates.")
    else:
        print(f"{counter} restraint centers are outside the sigma range of the starting coordinates.")
        print(output_string)

def set_default(sea_map, key, default):
    """
        Analagous to dict.setdefault.
    """
    sea_map[key] = sea_map[key] if key in sea_map else default
    return sea_map[key]

def add_restraints_to_stgs(stgs: sea.Map, restraints: dict, overwrite = False):
    """
        Add the restraints to the first assign_forcefield stage.
        We use the assign_forcefield stage so the terms will not be ignored
        or altered by `restraints.existing` parameters in later stages.

        Parameters:
        ----------
        stgs : sea.Map
            The stages of the MSJ as a sea.Map object.
        restraints : dict
            A dictionary containing the restraints to be added to the MSJs. The keys should be 
            the type of restraints (e.g. 'posre_fbhw', 'posre_harm') and the values should be 
            lists of dictionaries containing the restraint information. Each dictionary should 
            contain the following keys: 'atoms', 'ref', 'force_constants', and 'sigma' (optional).
        overwrite : bool, optional
            Whether to overwrite existing restraints in the force field assignment stage of the MSJs.
            Default is False.   

        Returns:
        -------
        None
    """
    for stg in stgs['stage']:
        # Make sure that the system is not recentered at the building stage
        # because restraint center coordinates are defined explicitly
        if stg.__NAME__ == "build_geometry":
            stg.rezero_system = sea.Atom(False)
        # Add the restraints to the forcefield stage
        if stg.__NAME__ == "assign_forcefield":
            restraint_map = set_default(stg, 'restraints', sea.Map())
            if overwrite:
                old_restraints = stg.restraints.existing = sea.Atom('ignore')
                new_restraints = stg.restraints.new = sea.List()
            else:
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
        # Disable the restraint correction
        if stg.__NAME__ == "fep_analysis":
            stg.correct_restr = sea.Atom(False)

def add_restraints_to_graph_msjs(graph_id, restraints, overwrite=False):
    """
        Add restraints to the forcefield stage of the MSJs in a GraphDB graph.

        Parameters:
        ----------
        graph_id : str
            The ID of the graph in GraphDB.
        restraints : dict
            A dictionary containing the restraints to be added to the MSJs. The keys should be 
            the type of restraints (e.g. 'posre_fbhw', 'posre_harm') and the values should be 
            lists of dictionaries containing the restraint information. Each dictionary should 
            contain the following keys: 'atoms', 'ref', 'force_constants', and 'sigma' (optional).
        overwrite : bool, optional
            Whether to overwrite existing restraints in the force field assignment stage of the MSJs.
            Default is False.

        Returns:
        -------
        None
    """
    graph = driver.Graph.load(graph_id)
    for edge in graph.edges():
        edge_update = {}
        for leg in edge.get_legs():
            if leg.leg_type != 'complex':
                continue
            msj = add_restraints_to_msj_string(leg.msj, restraints, overwrite=overwrite)
            edge_update.update({
                leg.leg_type: {
                    'msj': driver.to_string(
                        base64.b64encode(driver.to_bytes(msj)))
                }
            })
        edge.update(edge_update)

def add_restraints_to_msj_string(msj_content, restraints, overwrite=False) -> str:
    """
        Add restraints to the forcefield stage of a Desmond MSJ string. 

        Parameters:
        ----------
        msj_content : str
            The content of the MSJ as a string.
        restraints : dict
            A dictionary containing the restraints to be added to the MSJs. The keys should be 
            the type of restraints (e.g. 'posre_fbhw', 'posre_harm') and the values should be 
            lists of dictionaries containing the restraint information. Each dictionary should 
            contain the following keys: 'atoms', 'ref', 'force_constants', and 'sigma' (optional).
        overwrite : bool, optional
            Whether to overwrite existing restraints in the force field assignment stage of the MSJs.
            Default is False.

        Returns:
        -------
        str
            The updated MSJ content as a string.
    """
    stgs = cmj.msj2sea('', msj_content)
    add_restraints_to_stgs(stgs, restraints, overwrite=overwrite)
    stgs.add_tag("setbyuser")
    msj_content = cmj.write_sea2msj(stgs.stage, to_str=True)
    msj_parser = parse(string=msj_content)
    return str(msj_parser)

def submit_graphdb_job_with_restraints(fmp_or_pv_file, yaml_file, restraints, overwrite_restraints=False):
    """
        This function submits a job to GraphDB using the provided FMP or PV file and YAML file. 
        It adds the specified restraints to the first `assign_forcefield` stage of the MSJs 
        (Molecular Simulation Jobs) and optionally overwrites/replaces existing restraints.

        Parameters:
        ----------
        fmp_or_pv_file : str
            The path to the FMP or PV file that will be used to generate the graph.
        yaml_file : str
            The path to the YAML file that will be used to generate the graph.
        restraints : dict
            A dictionary containing the restraints to be added to the MSJs. The keys should be 
            the type of restraints (e.g. 'posre_fbhw', 'posre_harm') and the values should be 
            lists of dictionaries containing the restraint information. Each dictionary should 
            contain the following keys: 'atoms', 'ref', 'force_constants', and 'sigma' (optional).
        overwrite_restraints : bool, optional
            Whether to overwrite existing restraints in the force field assignment stage of the MSJs.
            Default is False.

        Returns:
        -------
        None
    """
    url = 'https://graphdb.schrodinger.com' # Make sure to use this "private"/advanced URL
    jws.login(web_services_addr=url)
    gid = jws.submit(fmp_or_pv_file, yaml_file, graph_generation_only=True)
    add_restraints_to_graph_msjs(gid, restraints, overwrite=overwrite_restraints)
    jws.requeue(gid)
