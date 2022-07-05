import sys, os
import pickle
import math
import re
import tarfile
import tempfile
import shutil
from operator import itemgetter

from schrodinger.application.prime.packages import executeparallel
from schrodinger.application.prime.packages import primeparser
from schrodinger import structure
from schrodinger.protein import residue
from schrodinger.structutils import analyze, rmsd
import schrodinger.application.desmond.cms as cms
from schrodinger.application.desmond.packages import topo, traj_util

import numpy as np
import pylab
from matplotlib import patches
from matplotlib.collections import PatchCollection

usage = ""
version = ""

class TrajectoryRMSF(executeparallel.ExecuteParallel):
    def __init__(self, script_file):
        executeparallel.ExecuteParallel.__init__(self, script_file,
            usage, [], version=version)

    def addParserArguments(self, parser):
        parser.add_argument("trj",
                action=primeparser.StoreFile,
                help="output CMS")
        parser.add_argument("reference",
                action=primeparser.StoreFile,
                help="Reference structure to align to")
        parser.add_argument("--pickle",
                action=primeparser.StoreFile,
                help="Pickle file containing data, only plot")
        parser.add_argument("--asl", action=primeparser.StoreAsl,
                default="protein AND backbone",
                help='ASL for loop')
        parser.add_argument("--fitasl", action=primeparser.StoreAsl,
                default="protein AND backbone",
                help='ASL for fitting')
        parser.add_argument("--dmorph", choices=('a', 'b'),
                default=None,
                help="Side to use in a deactivated morphing trajectory")
        parser.add_argument('--min', type=float,
                default=0, help='Minimum RMSF for plot')
        parser.add_argument('--max', type=float,
                default=None, help='Maximum RMSF for plot')
        parser.add_argument("--structure",
                action='store_true', default=False,
                help="Output average structure")
        parser.add_argument("--use_ref",
                action='store_true',
                default=False,
                help="Use the reference structure instead of the average structure.")

    def runBackend(self, opts, istage):
        if opts.pickle:
            with open(opts.pickle, 'rb') as fh:
                data = pickle.load(fh)
            extracted_pos = data['extracted_pos']
            rmsf_reference_pos = data['rmsf_reference_pos']
            residue_rmsf = data['residue_rmsf']
            residues = data['residues']
            fit_rmsd = data['fit_rmsd']
        else:
            residue_rmsf, residues, extracted_pos, rmsf_reference_pos, fit_rmsd = self._parseTarball(opts)
            # Save
            output_fname = opts.jobname + '.pickle'
            with open(output_fname, 'wb') as fh:
                data = {'extracted_pos': extracted_pos,
                        'rmsf_reference_pos': rmsf_reference_pos,
                        'residue_rmsf': residue_rmsf,
                        'residues': residues,
                        'fit_rmsd': fit_rmsd}
                pickle.dump(data, fh)
            self.addOutputFile(output_fname)
        self._plotRMSF(opts, residue_rmsf, residues)
        self._plotRMSD(opts, fit_rmsd)
        if opts.structure:
            self._outputStructures(opts, rmsf_reference_pos, residues, residue_rmsf)

    def _parseTarball(self, opts):
        '''
        '''
        # Find a good reference structure
        with structure.StructureReader(self.getPath(opts.reference)) as structures:
            for st_idx, ref_st in enumerate(structures):
                ref_fitasl_alist = analyze.evaluate_asl(ref_st, opts.fitasl)
                if len(ref_fitasl_alist) == 0:
                    raise UserWarning('{0:s} structure {1:d} has no matching atoms for {2:s}'.format(opts.reference, st_idx, opts.fitasl))
                print('"{0:s}" matches {1:d} atoms in {2:s} structure {3:d}'.format(opts.fitasl, len(ref_fitasl_alist), opts.reference, st_idx))
                ref_asl_alist = analyze.evaluate_asl(ref_st, opts.asl)
                if len(ref_asl_alist) == 0:
                    raise UserWarning('{0:s} structure {1:d} has no matching atoms for {2:s}'.format(opts.reference, st_idx, opts.asl))
                # Break down into residues
                residues = tuple(sorted({ref_st.atom[aid].resnum for aid in ref_asl_alist}))
                print(f'"{opts.asl:s}" matches {len(ref_asl_alist):d} atoms ({len(residues)} residues) in {opts.reference:s} structure {st_idx:d}')
                break
        #######################################
        # Data structure: {R} F x A x 3       #
        # R: res1, res2 ... resR              #
        # F: frame                            #
        # A: atom1 atom2 ... atomA            #
        # 3: posX posY posZ                   #
        #######################################
        # Don't know the number of frames yet
        fit_rmsd = []
        extracted_pos = {}
        fit_idx = None
        # Process each replica
        _, cms_model, trj = traj_util.read_cms_and_traj(self.getPath(opts.trj))
        fsys_ct = cms_model.fsys_ct.copy()
        # Fitting
        resnum_alist_0idx = {}
        ref_resnum_alist_0idx = {}
        trj_fitasl = None
        trj_asl = None
        if opts.dmorph:
            if opts.dmorph == 'a':
                trj_fitasl = f"{opts.fitasl} AND atom.i_fep_side 1"
                trj_asl = f"{opts.asl} AND atom.i_fep_side 1"
            elif opts.dmorph == 'b':
                trj_fitasl = f"{opts.fitasl} AND atom.i_fep_side 2"
                trj_asl = f"{opts.asl} AND atom.i_fep_side 2"
            else:
                raise UserWarning(f"invalid dmorph option {opts.dmorph}")
        else:
            trj_fitasl = opts.fitasl
            trj_asl = opts.asl
        trj_fitasl_alist = analyze.evaluate_asl(cms_model, trj_fitasl)
        if len(trj_fitasl_alist) == 0:
            raise UserWarning('trajectory has no matching atoms for {0:s}'.format(trj_fitasl))
        elif len(trj_fitasl_alist) != len(ref_fitasl_alist):
            raise UserWarning('{0:d} trj atoms vs {1:d} ref atoms for "{2:s}"'.format(len(trj_fitasl_alist), len(ref_fitasl_alist), opts.fitasl))
        print('"{0:s}" matches {1:d} atoms in {2:s}'.format(trj_fitasl, len(trj_fitasl_alist), opts.trj))
        trj_asl_alist = analyze.evaluate_asl(cms_model, trj_asl)
        if len(trj_asl_alist) == 0:
            raise UserWarning('trajectory has no matching atoms for {0:s}'.format(trj_asl))
        elif len(trj_asl_alist) != len(ref_asl_alist):
            raise UserWarning('{0:d} trj atoms vs {1:d} ref atoms for "{2:s}"'.format(len(trj_asl_alist), len(ref_asl_alist), opts.asl))
        print('"{0:s}" matches {1:d} atoms in {2:s}'.format(trj_asl, len(trj_asl_alist), opts.trj))
        # Setup the data
        for resnum in residues:
            trj_resnum_asl_alist = analyze.evaluate_asl(cms_model, f"{trj_asl} AND res.num {resnum}")
            # The alists are 1-indexed, but frame_pos needs to be 0-indexed
            trj_resnum_alist_0idx = [val-1 for val in trj_resnum_asl_alist]
            if resnum not in extracted_pos:
                extracted_pos[resnum] = np.zeros((len(trj), len(trj_resnum_alist_0idx), 3))
            resnum_alist_0idx[resnum] = trj_resnum_alist_0idx
            # And for the reference structure, if necessary
            if opts.use_ref:
                ref_resnum_asl_alist = analyze.evaluate_asl(ref_st, f"{trj_asl} AND res.num {resnum}")
                # The alists are 1-indexed, but frame_pos needs to be 0-indexed
                ref_resnum_alist_0idx[resnum] = [val-1 for val in ref_resnum_asl_alist]
        # Compute the RMSD
        for iframe, frame in enumerate(trj):
            # Update the ct positions to the traj positions
            # We have to do this because traj includes pseudo-atoms
            topo.update_ct(fsys_ct, cms_model, frame)
            # First fit the entire CT since we have at least two molecules
            rms_fit = rmsd.superimpose(ref_st, ref_fitasl_alist,
                    fsys_ct, trj_fitasl_alist, move_which=rmsd.CT)
            fit_rmsd.append(rms_fit)
            # Then extract the positions, assuming the same ordering as in the reference st
            frame_pos = fsys_ct.getXYZ(copy=True)
            for resnum, alist_0idx in resnum_alist_0idx.items():
                resnum_pos = extracted_pos[resnum]
                resnum_pos[iframe, :] = frame_pos[alist_0idx, :]
        # Compute the RMSF reference (usually mean)
        #######################################
        # Data structure: {R} 1 x A           #
        # R: res1, res2, res3 ... resR        #
        # 1: placeholder for frame            #
        # A: atom1, atom2 ... atomA (per res) #
        #######################################
        rmsf_reference_pos = {}
        if opts.use_ref:
            for res_num in residues:
                num_atoms = len(ref_resnum_alist_0idx[res_num])
                ref_pos = ref_st.getXYZ()[ref_resnum_alist_0idx[res_num]]
                rmsf_reference_pos[res_num] = ref_pos
        else:
            for res_num, res_pos in extracted_pos.items():
                avg_shape = list(res_pos.shape)
                avg_shape[0] = 1
                avg_pos = res_pos.mean(axis=0).reshape(avg_shape)
                rmsf_reference_pos[res_num] = avg_pos
        # Compute the actual RMSFs
        #######################################
        # Data structure: R                   #
        # R: res1, res2, res3 ... resR        #
        #######################################
        residue_rmsf = np.zeros(len(residues))
        for res_idx, (res_num, res_pos) in enumerate(extracted_pos.items()):
            ref_pos = rmsf_reference_pos[res_num]
            offset = np.square(res_pos - ref_pos).sum(axis=2)
            rmsf = np.sqrt(offset.mean(axis=(0,1))) # Average over frames AND atoms
            residue_rmsf[res_idx] = rmsf
        # And make RMSD a numpy array
        #######################################
        # Data structure: F                   #
        # R: frame                            #
        #######################################
        fit_rmsd = np.array(fit_rmsd)
        return residue_rmsf, residues, extracted_pos, rmsf_reference_pos, fit_rmsd

    def _plotRMSF(self, opts, residue_rmsf, residues):
        # Find a good reference structure
        with structure.StructureReader(self.getPath(opts.reference)) as structures:
            for st_idx, st in enumerate(structures):
                fitasl_alist = analyze.evaluate_asl(st, opts.fitasl)
                if len(fitasl_alist) == 0:
                    raise UserWarning('{0:s} structure {1:d} has no matching atoms for {2:s}'.format(opts.reference, st_idx, opts.fitasl))
                print('"{0:s}" matches {1:d} atoms in {2:s} structure {3:d}'.format(opts.fitasl, len(fitasl_alist), opts.reference, st_idx))
                asl_alist = analyze.evaluate_asl(st, opts.asl)
                if len(asl_alist) == 0:
                    raise UserWarning('{0:s} structure {1:d} has no matching atoms for {2:s}'.format(opts.reference, st_idx, opts.asl))
                # Break down into residues
                current_residues = tuple(sorted({st.atom[aid].resnum for aid in asl_alist}))
                print(f'"{opts.asl:s}" matches {len(asl_alist):d} atoms ({len(residues)} residues) in {opts.reference:s} structure {st_idx:d}')
                if current_residues != residues:
                    raise UserWarning("reference residues did not match extracted residues")
                break
            else:
                raise UserWarning("No appropriate structure was found")
        # Extract the secondary structure
        secondary = structure_to_secondary_structure(st)
        ss_colors = {
                'Helix': 'b',
                'Strand': 'g'
                }
        # Plot
        fig = pylab.figure()
        ax = pylab.gca()
        pylab.plot(residues, residue_rmsf)
        add_secondary_structure_to_axis(st, pylab.gca(), residues[0], residues[-1])
        # Print out the span of the >3A RMSF
        high_indices = np.where(residue_rmsf > 3)[0]
        breaks = np.nonzero(np.not_equal(high_indices[1:], high_indices[:-1] + 1))[0]
        start_idx = 0
        for end_idx in breaks:
            print(f"Residues {residues[high_indices[start_idx]]} to {residues[high_indices[end_idx]]} had >3.0A RMSF")
            pylab.axvspan(residues[high_indices[start_idx]], residues[high_indices[end_idx]], color='g', alpha=0.25)
            start_idx = end_idx+1
        #print(f"Span {residues[high_indices[start_idx]]} to {residues[high_indices[-1]]} had >3.0A RMSF")
        #pylab.axvspan(residues[high_indices[start_idx]], residues[high_indices[-1]], color='g', alpha=0.25)
        pylab.xlabel('Residue')
        pylab.ylabel('RMSF (A)')
        ax.set_ylim(bottom=opts.min)
        if opts.max:
            ax.set_ylim(top=opts.max)
        fig_fname = f"{opts.jobname:s}_rmsf.png"
        pylab.savefig(fig_fname)
        self.addOutputFile(fig_fname)
        pylab.close()
    
    def _plotRMSD(self, opts, fit_rmsd):
        pylab.figure()
        ax = pylab.gca()
        pylab.plot(fit_rmsd)
        pylab.xlabel('Frame')
        pylab.ylabel('RMSD (A)')
        ax.set_ylim(bottom=opts.min)
        if opts.max:
            ax.set_ylim(top=opts.max)
        pylab.title(f"Fit RMSD")
        fig_fname = f"{opts.jobname:s}_fit_rmsd.png"
        pylab.savefig(fig_fname)
        self.addOutputFile(fig_fname)
        pylab.close()

    def _outputStructures(self, opts, rmsf_reference_pos, residues, residue_rmsf):
        '''Output the average positions superimposed on the reference structure.

        TODO RMSF stored in the beta factors.
        '''
        # Find a good reference structure
        with structure.StructureReader(self.getPath(opts.reference)) as structures:
            for st_idx, ref_st in enumerate(structures):
                ref_fitasl_alist = analyze.evaluate_asl(ref_st, opts.fitasl)
                ref_asl_alist = analyze.evaluate_asl(ref_st, opts.asl)
                break
        # Remove the atoms on ASL residues that are not covered by the ASL (will look super funky otherwise)
        residue_asl_alist = analyze.evaluate_asl(ref_st, f"fillres ({opts.asl})")
        delete_alist = set(residue_asl_alist).difference(ref_asl_alist)
        ref_st.deleteAtoms(delete_alist)
        # Set beta factors to 0 for everything
        for atom in ref_st.atom:
            atom.temperature_factor = 0.0
        # Update positions for each replica
        resnum_alists = {}
        for resnum in residues:
            #FIXME: this assumes ref_st and the trajectory return the same order of atoms for each residue
            trj_resnum_asl_alist = analyze.evaluate_asl(ref_st, f"{opts.asl} AND res.num {resnum}")
            resnum_alists[resnum] = trj_resnum_asl_alist
        num_replicas = next(iter(rmsf_reference_pos.values())).shape[0]
        for replica_i in range(num_replicas):
            for res_idx, (resnum, ref_pos) in enumerate(rmsf_reference_pos.items()):
                # Grab the right replica so we have A x 3
                atom_xyzs = ref_pos[replica_i, 0, :, :]
                for atom_idx, atom_xyz in zip(resnum_alists[resnum], atom_xyzs):
                    ref_st.atom[atom_idx].xyz = atom_xyz
                    # Temperature factor in Maestro goes from Blue (<5) to Red (>67)
                    # So we scale RMSF by a factor of 20
                    ref_st.atom[atom_idx].temperature_factor = residue_rmsf[replica_i, res_idx] * 20
            # Write out structure
            fname = f"{opts.jobname}_ref_pos_replica{replica_i:03d}.mae"
            ref_st.write(fname)
            self.addOutputFile(fname)

    def _plotReplicaReferenceRMSF(self, opts, rmsf_reference_pos, residues):
        # Find a good reference structure
        with structure.StructureReader(self.getPath(opts.reference)) as structures:
            for st_idx, ref_st in enumerate(structures):
                ref_fitasl_alist = analyze.evaluate_asl(ref_st, opts.fitasl)
                if len(ref_fitasl_alist) == 0:
                    raise UserWarning('{0:s} structure {1:d} has no matching atoms for {2:s}'.format(opts.reference, st_idx, opts.fitasl))
                print('"{0:s}" matches {1:d} atoms in {2:s} structure {3:d}'.format(opts.fitasl, len(ref_fitasl_alist), opts.reference, st_idx))
                ref_asl_alist = analyze.evaluate_asl(ref_st, opts.asl)
                if len(ref_asl_alist) == 0:
                    raise UserWarning('{0:s} structure {1:d} has no matching atoms for {2:s}'.format(opts.reference, st_idx, opts.asl))
                # Break down into residues
                ref_residues = tuple(sorted({ref_st.atom[aid].resnum for aid in ref_asl_alist}))
                print(f'"{opts.asl:s}" matches {len(ref_asl_alist):d} atoms ({len(residues)} residues) in {opts.reference:s} structure {st_idx:d}')
                if residues != ref_residues:
                    raise UserWarning("reference residues did not match extracted residues")
                break
        # Extract the secondary structure
        secondary = structure_to_secondary_structure(ref_st)
        ss_colors = {
                'Helix': 'b',
                'Strand': 'g'
                }
        # rmsf_reference_pos
        #######################################
        # Data structure: {R} N x 1 x A       #
        # R: res1, res2, res3 ... resR        #
        # N: replica1, replica2 ... replicaN  #
        # 1: placeholder for frame            #
        # A: atom1, atom2 ... atomA (per res) #
        #######################################
        num_residues = len(rmsf_reference_pos)
        num_replicas = next(iter(rmsf_reference_pos.values())).shape[0]
        # Compute the RMSF
        #######################################
        # Data structure: R                   #
        # R: res1, res2, res3 ... resR        #
        #######################################
        reference_pos_rmsf_vs_first = np.zeros(num_residues)
        reference_pos_rmsf_vs_last = np.zeros(num_residues)
        for res_idx, (res_num, res_pos) in enumerate(rmsf_reference_pos.items()):
            # RMSF vs first average structure
            offset = np.square(res_pos[1:, :] - res_pos[0, :]).sum(axis=3)
            rmsf = np.sqrt(offset.mean(axis=(0,1,2))) # Average over replicas, frames AND atoms
            reference_pos_rmsf_vs_first[res_idx] = rmsf
            # RMSF vs last average structure
            offset = np.square(res_pos[:-1, :] - res_pos[-1, :]).sum(axis=3)
            rmsf = np.sqrt(offset.mean(axis=(0,1,2))) # Average over replicas, frames AND atoms
            reference_pos_rmsf_vs_last[res_idx] = rmsf
        for label, rmsf in (('first', reference_pos_rmsf_vs_first),
                            ('last', reference_pos_rmsf_vs_last)):
            fig = pylab.figure()
            ax = pylab.gca()
            pylab.plot(residues, rmsf)
            add_secondary_structure_to_axis(ref_st, pylab.gca(), residues[0], residues[-1])
            pylab.xlabel('Residue')
            pylab.ylabel('RMSF (A)')
            ax.set_ylim(bottom=opts.min)
            if opts.max:
                ax.set_ylim(top=opts.max)
            fig_fname = f"{opts.jobname:s}_reference_rmsf_vs_{label}.png"
            pylab.savefig(fig_fname)
            self.addOutputFile(fig_fname)
            pylab.close()


def structure_to_secondary_structure(st):
    '''Returns a dictionary of SS_TYPE: [(res_start, res_end)...]'''
    secondary = {}
    prev_ss = None
    start_resnum = None
    prev_resnum = None
    for res in st.residue:
        if res.secondary_structure != prev_ss:
            if start_resnum is not None:
                ss_runs = secondary.setdefault(residue.SSA_TT_MAP[prev_ss], [])
                ss_runs.append((start_resnum, prev_resnum))
            prev_ss = res.secondary_structure
            start_resnum = res.resnum
        prev_resnum = res.resnum
    return secondary

def add_secondary_structure_to_axis(st, ax, first_resid, last_resid):
    # Extract the secondary structure
    xlim = ax.get_xlim()
    secondary = structure_to_secondary_structure(st)
    CENTER_Y = 0.2
    HELIX_HEIGHT = 0.15
    STRAND_HEIGHT = 0.1
    CHAIN_HEIGHT = 0.05
    to_add = []
    to_add.append(patches.Rectangle((first_resid, CENTER_Y-CHAIN_HEIGHT/2), (last_resid-first_resid), CHAIN_HEIGHT, color='k'))
    # Add the secondary structure shapes
    for start_resnum, end_resnum in secondary.get('Helix', ()):
        patch = patches.Rectangle((start_resnum, CENTER_Y-HELIX_HEIGHT/2), (end_resnum-start_resnum), HELIX_HEIGHT, color='y')
        to_add.append(patch)
    for start_resnum, end_resnum in secondary.get('Strand', ()):
        patch = patches.Arrow(start_resnum, CENTER_Y, (end_resnum-start_resnum), 0, width=STRAND_HEIGHT, color='g')
        to_add.append(patch)
    collection = PatchCollection(to_add)
    ax.add_collection(collection)
    ax.set_xlim(xlim)

if __name__ == "__main__":
    TrajectoryRMSF(sys.argv[0]).launch(sys.argv[1:])
