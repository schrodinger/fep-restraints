import os
import shutil
from pathlib import Path
import schrodinger
import schrodinger.structure as structure
from schrodinger.application.desmond import constants
from schrodinger.application.desmond.packages import cui
#from schrodinger.application.desmond.starter.generator.abfep import prepare_inputs, generate
from schrodinger.application.desmond.starter.generator.fep_absolute_binding import prepare_inputs, generate
from schrodinger.forcefield.custom_params import create_archive_from_oplsdir, merge_oplsdir, upgrade_oplsdir


def write_abfep_restraints_job_script(
    job_name, rest_file, align_ref, atom_sel, fep_force_const, md_force_const, scaling, 
    fep_sim_time, md_sim_time, host, subhost, maxjob, retries, 
    salt=None, membrane=False, membrane_type='POPC',
    ff='OPLS4', opls=None, ffhost=None, bottom_width=None, 
    project=None, account=None, qarg=None, seed=2007
    ):

    # Get the directory where the code is
    code_file_path = os.path.realpath(__file__)
    code_directory = os.path.dirname(code_file_path)
    
    with open(job_name+'.sh', mode = 'wt', encoding = 'utf-8') as f:
        # Preparation step
        f.write('# Run the preparation step\n'
                f'$SCHRODINGER/fep_absolute_binding \\\n  "{job_name}.fmp" \\\n'
                f'  -JOBNAME "{job_name}" \\\n  -prepare -ppj 4 -maxjob {maxjob} -ensemble muVT \\\n'
                f'  -fep-sim-time {fep_sim_time} -md-sim-time {md_sim_time} -ff {ff} -seed {seed} \\\n'
                f'  -HOST "{host}" \\\n  -SUBHOST "{subhost}"')
        if salt is not None:
            f.write(f' \\\n  -salt {salt} ')
        if membrane:
            f.write(f' \\\n  -membrane \\\n  -membrane-type {membrane_type}')
        if opls is not None: 
            f.write(f' \\\n  -OPLSDIR "{opls}" ')
        if ffhost is not None:
            f.write(f' \\\n  -ffbuilder \\\n  -ff-host "{ffhost}"')
        f.write('\n')
        # Restraints for FEP of complex
        f.write('\n# Add restraints to the FEP simulations of the complex\n'
                'for MSJ in *_complex.msj; do\n  OLD="${MSJ}.old"\n  mv "$MSJ" "$OLD"\n'
                f'  $SCHRODINGER/run {code_directory}/add_sigma_restraints.py \\\n'
                f'    "{rest_file}" "$OLD" "$MSJ" -f {fep_force_const} -s {scaling} \\\n'
                f'    -r "{align_ref}" \\\n    -a "{atom_sel}"')
        if bottom_width is not None: 
            f.write(f' \\\n    --bw {bottom_width} ')
        f.write('\ndone\n')
        # Restraints for initial MD
        f.write('\n# Add restraints to the initial MD run.\n'
                'for MSJ in *_md.msj; do\n  OLD="${MSJ}.old"\n  mv "$MSJ" "$OLD"\n'
                f'  $SCHRODINGER/run {code_directory}/add_sigma_restraints.py \\\n'
                f'    "{rest_file}" "$OLD" "$MSJ" -f {md_force_const} -s {scaling} \\\n'
                f'    -r "{align_ref}" \\\n    -a "{atom_sel}"')
        if bottom_width is not None: 
            f.write(f' \\\n    --bw {bottom_width} ')
        f.write('\ndone\n')
	# Frequency for MD trajectory
        f.write('\nsed -i \'s/interval = 3.6/interval = 36.0/g\' *_md.msj\n')
        # Submission
        f.write('\n# Submit the job to the cluster.\n'
                f'$SCHRODINGER/utilities/multisim \\\n  "{job_name}_pv.maegz" \\\n'
                f'  -JOBNAME "{job_name}" \\\n  -m "{job_name}.msj" \\\n  -o "{job_name}-out.mae" \\\n'
                f'  -HOST "{host}" -SUBHOST "{subhost}" \\\n'
                f'  -maxjob {maxjob} -RETRIES {retries}')
        # Construct the QARG string (including project or account if given)
        qarg_out = qarg
        if project is not None:
            if qarg_out is None:
                qarg_out = f'-P {project}'
            else:
                qarg_out += f' -P {project}'
        if account is not None:
            if qarg_out is None:
                qarg_out = f'-A {account}'
            else:
                qarg_out += f' -A {account}'
        if qarg_out is not None:
            f.write(f' \\\n  -QARG "{qarg_out}"')
        if opls is not None: 
            f.write(f' \\\n  -OPLSDIR "{opls}"')
        f.write('\n')


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Extract individual pv files from one FEP pv file.")
    parser.add_argument('mae_file', help="File with input structures for FEP")
    parser.add_argument('-j','--job-name', help="Name for the output directory")
    parser.add_argument('-r','--rest-file', type=str, default=None, help='File with the restraint center definitions and optionally restraint bottom widths.')
    parser.add_argument('-s','--scaling-factor', type=float, default=1.0, help='Scaling factor for the bottom width of the restraints. Default: 1.0')
    parser.add_argument('-b','--bottom-width', type=float, default=None, help='If the restraint file does not contain bottom widths, this value is used for all restraints. If neither this nor a bottom width in the restraint file is given, a default of 0 is used.')
    parser.add_argument('--align-ref', type=str, default=None, help='Alignment reference file. If not given, the input file is used as reference.')
    parser.add_argument('--asl', dest='atom_sel', type=cui.validate_and_enclose_asl, default='all', help='Atom selection for the restrained atoms. Default: all.')
    parser.add_argument('--md-force-const', type=float, default=0.1, help='Force constant for the restraints in the initial MD stage. Default: 0.1')
    parser.add_argument('--fep-force-const', type=float, default=1.0, help='Force constant for the restraints in the FEP stage. Default: 1.0')
    parser.add_argument('--host', type=str, default='driver', help='Host for the job submission. Default: driver')
    parser.add_argument('--subhost', type=str, default='gpu', help='Subhost for the job submission. Default: gpu')
    parser.add_argument('--project', type=str, default=None, help='Project name for the job submission (only on Bolt).')
    parser.add_argument('--account', type=str, default=None, help='Account name for the job submission (only on Ada).')
    parser.add_argument('--qarg', type=str, default=None, help='Additional arguments for the job submission.')
    parser.add_argument('--maxjob', type=int, default=0, help='Maximum number of jobs to run in parallel. Default: 0 (no limit).')
    parser.add_argument('--retries', type=int, default=5, help='Number of retries for failed jobs. Default: 5.')
    parser.add_argument('--ff', choices=['OPLS4', 'OPLS5'], default='OPLS4', help='Specify the force field to use. Default: OPLS4.')
    parser.add_argument('--oplsdir', type=str, default=None, help='Path to the OPLS force field directory. Default: None.')
    parser.add_argument('--ffbuilder', action='store_true', help='Start a force field builder job for missing parameters.')
    parser.add_argument('--ffhost', type=str, default='cpu', help='Host for the force field builder job. Default: cpu.')
    parser.add_argument('--md-sim-time', type=float, default=2000, help='MD simulation time in ps. Default: 2000.')
    parser.add_argument('--fep-sim-time', type=float, default=10000, help='FEP simulation time in ps. Default: 10000.')
    parser.add_argument('--salt', type=float, default=None, help='Salt concentration in Molar. Default: None (no salt).')
    parser.add_argument('--membrane', action='store_true', help='Indicates the model system is a membrane-bound protein system, such as the GPCR or an ion channel. If membrane/water components are not provided, the bilayer specified by --membrane-type will be added and equilibrated. The protein coordinates should be OPM-compatible. If the membrane is included in the inputs, then the equilibration will be skipped. Default: False.')
    parser.add_argument('--membrane-type', choices=['POPC', 'DPPC', 'DMPC', 'POPE'], default='POPC', help='Specify the type of membrane to use. Default: POPC.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite the job directory if it already exists.')
    parser.add_argument('--lig-rest', action='store_true', help='Use the ligand restraints from the input file. Default: False.')
    parser.add_argument('--seed', type=int, default=2007, help='Random seed for the job. Default: 2007.')
    args = parser.parse_args()

    # Remove FFB host if no ff is generated
    if not args.ffbuilder:
        ffhost = None
    # File paths from input 
    align_ref = args.align_ref
    if args.align_ref is None:
        align_ref = args.mae_file
    mae_path = Path(args.mae_file)
    job_path = Path(args.job_name)
    rest_path = Path(args.rest_file)
    align_path = Path(align_ref)
    # Create the directory
    if os.path.isdir(job_path):
        print(f'Directory {job_path} already exists!')
        if args.overwrite:
            print(f'It will be overwritten.')
            shutil.rmtree(job_path)
        else:
            exit
    os.makedirs(job_path)
    # Copy the start file into it
    pvfile = f'{args.job_name}-in_pv{mae_path.suffix}' 
    shutil.copy(mae_path, Path(job_path, pvfile))
    # Copy the restraints file into it
    rsfile = f'{args.job_name}-restr.cms'
    shutil.copy(rest_path, Path(job_path, rsfile))
    # Copy the restraints alignment file into it
    alfile = f'{args.job_name}-align{align_path.suffix}' 
    shutil.copy(align_path, Path(job_path, alfile))
    # Copy the (optional) force field dir into it
    if args.oplsdir is not None:
        shutil.copy(Path(args.oplsdir), job_path)
        upgrade_oplsdir(args.job_name)
        release = schrodinger.get_release_name().replace('-','_')
        opls = f'custom_{release}.opls'
    else:
        opls = None

    # Enter the job directory
    os.chdir(job_path)
    # Create the fmp file
    _, _graph = prepare_inputs(Path(pvfile), args.lig_rest, Path(args.job_name))
    # Create the job submission script
    write_abfep_restraints_job_script(
        args.job_name, rsfile, alfile, args.atom_sel, 
        args.fep_force_const, args.md_force_const, args.scaling_factor,
        args.fep_sim_time, args.md_sim_time, args.host, args.subhost, 
        args.maxjob, args.retries, 
        salt=args.salt,
        membrane=args.membrane,
        membrane_type=args.membrane_type,
        ff=args.ff,
        opls=opls, 
        ffhost=args.ffhost, 
        bottom_width=args.bottom_width,
        project=args.project,
        account=args.account,
        qarg=args.qarg,
        seed=args.seed
    )        
    # Leave the job directory
    os.chdir('..')

