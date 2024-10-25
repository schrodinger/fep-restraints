
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
    '''Returns an ASL that (should) uniquely identify the provided atom.
    '''
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
    '''Take a structure and create position restraint (posre) terms from it.

    If the matched atom has a r_desmond_sigma property it will be assigned to posre_fbhw, otherwise posre_harm.
    The coordinates of the atom are taken as the center of the restraint.
    '''
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
    '''Read the sigma values from a Desmond structure file and create restraints.'''
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

def add_restraints_to_graph_msjs(graph_id, restraints):
    '''Add restraints to the MSJs of a graph.'''
    graph = driver.Graph.load(graph_id)
    for edge in graph.edges():
        if edge.subjob not in ["fep", ""]:
            continue
        edge_update = {}
        for leg in edge.get_legs():
            msj = add_restraints_to_msj_string(leg.msj, restraints)
            edge_update.update({
                leg.leg_type: {
                    'msj': driver.to_string(
                        base64.b64encode(driver.to_bytes(msj)))
                }
            })
        edge.update(edge_update)

def add_restraints_to_msj_string(msj_content, restraints) -> str:
    '''Add restraints to the forcefield stage of a Desmond MSJ string.'''
    stgs = cmj.msj2sea('', msj_content)
    add_restraints_to_stgs(stgs, restraints)
    stgs.add_tag("setbyuser")
    msj_content = cmj.write_sea2msj(stgs.stage, to_str=True)
    msj_parser = parse(string=msj_content)
    return str(msj_parser)

def submit_graphdb_job_with_restraints(fmp_or_pv_file, yaml_file, restraints):
    url = 'https://graphdb.schrodinger.com' # Make sure to use this "private"/advanced URL
    jws.login(web_services_addr=url)
    gid = jws.submit(fmp_or_pv_file, yaml_file, graph_generation_only=True)
    add_restraints_to_graph_msjs(gid, restraints)
    jws.requeue(gid)
