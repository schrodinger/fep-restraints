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
    fep_sim_time, md_sim_time, host, subhost, project, maxjob, retries, 
    salt=None, opls=None, ffhost=None, bottom_width=None
    ):

    # Get the directory where the code is
    code_file_path = os.path.realpath(__file__)
    code_directory = os.path.dirname(code_file_path)
    
    with open(job_name+'.sh', mode = 'wt', encoding = 'utf-8') as f:
        # Preparation step
        f.write('# Run the preparation step\n'
                f'$SCHRODINGER/fep_absolute_binding \\\n  "{job_name}.fmp" \\\n'
                f'  -JOBNAME "{job_name}" \\\n  -prepare -ppj 4 -maxjob {maxjob} -ensemble muVT \\\n'
                f'  -fep-sim-time {fep_sim_time} -md-sim-time {md_sim_time} -seed 2007 \\\n'
                f'  -HOST "{host}" \\\n  -SUBHOST "{subhost}"')
        if salt is not None:
            f.write(f' \\\n  -salt {salt} ')
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
        # Submission
        f.write('\n# Submit the job to the cluster.\n'
                f'$SCHRODINGER/utilities/multisim \\\n  "{job_name}_pv.maegz" \\\n'
                f'  -JOBNAME "{job_name}" \\\n  -m "{job_name}.msj" \\\n  -o "{job_name}-out.mae" \\\n'
                f'  -HOST "{host}" -SUBHOST "{subhost}" \\\n'
                f'  -maxjob {maxjob} -RETRIES {retries}')
        if project is not None:
            f.write(f'\\\n  -QARG \"-P {project}\"')
        if opls is not None: 
            f.write(f'\\\n  -OPLSDIR "{opls}"')
        f.write('\n')


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Extract individual pv files from one FEP pv file.")
    parser.add_argument('mae_file', help="File with input structures for FEP")
    parser.add_argument('-j','--job-name', help="Name for the output directory")
    parser.add_argument('-r','--rest-file', type=str, default=None)
    parser.add_argument('--scaling-factor', type=float, default=1.0)
    parser.add_argument('--bottom-width', type=float, default=None)
    parser.add_argument('--align-ref', type=str, default=None)
    parser.add_argument('--asl', dest='atom_sel', type=cui.validate_and_enclose_asl, default='all')
    parser.add_argument('--md-force-const', type=float, default=0.1)
    parser.add_argument('--fep-force-const', type=float, default=1.0)
    parser.add_argument('--host', type=str, default='driver-4core-standard')
    parser.add_argument('--subhost', type=str, default='gpu-t4-4x-ondemand')
    parser.add_argument('--project', type=str, default=None)
    parser.add_argument('--maxjob', type=int, default=0)
    parser.add_argument('--retries', type=int, default=5)
    parser.add_argument('--oplsdir', type=str, default=None)
    parser.add_argument('--ffbuilder', action='store_true')
    parser.add_argument('--ffhost', type=str, default='compute-16core-64gb-ondemand')
    parser.add_argument('--md-sim-time', type=float, default=2000)
    parser.add_argument('--fep-sim-time', type=float, default=10000)
    parser.add_argument('--salt', type=float, default=None)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--lig-rest', action='store_true')
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
        opls = f'custom_{release}'
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
        args.project, args.maxjob, args.retries, 
        salt=args.salt,
        opls=opls, 
        ffhost=args.ffhost, 
        bottom_width=args.bottom_width
    )        
    # Leave the job directory
    os.chdir('..')

