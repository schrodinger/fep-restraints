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


def write_abfep_job_script(
    job_name, fep_sim_time, md_sim_time, host, subhost, project, maxjob, retries, 
    salt=None, ff='OPLS4', opls=None, ffhost=None, seed=2007
    ):

    # Get the directory where the code is
    code_file_path = os.path.realpath(__file__)
    code_directory = os.path.dirname(code_file_path)
    
    with open(job_name+'.sh', mode = 'wt', encoding = 'utf-8') as f:
        # Preparation step
        f.write('# Run the preparation step\n'
                f'$SCHRODINGER/fep_absolute_binding \\\n  "{job_name}.fmp" \\\n'
                f'  -JOBNAME "{job_name}" \\\n  -ppj 4 -maxjob {maxjob} -ensemble muVT \\\n'
                f'  -fep-sim-time {fep_sim_time} -md-sim-time {md_sim_time} -ff {ff} -seed {seed} \\\n'
                f'  -HOST "{host}" \\\n  -SUBHOST "{subhost}"')
        if salt is not None:
            f.write(f' \\\n  -salt {salt} ')
        if opls is not None: 
            f.write(f' \\\n  -OPLSDIR "{opls}" ')
        if ffhost is not None:
            f.write(f' \\\n  -ffbuilder \\\n  -ff-host "{ffhost}"')
        if project is not None:
            f.write(f' \\\n  -QARG \"-P {project}\"')
        f.write('\n')


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Extract individual pv files from one FEP pv file.")
    parser.add_argument('mae_file', help="File with input structures for FEP")
    parser.add_argument('-j','--job-name', help="Name for the output directory")
    parser.add_argument('--host', type=str, default='driver-4core-standard')
    parser.add_argument('--subhost', type=str, default='gpu-t4-4x-ondemand')
    parser.add_argument('--project', type=str, default=None)
    parser.add_argument('--maxjob', type=int, default=0)
    parser.add_argument('--retries', type=int, default=5)
    parser.add_argument('--ff', choices=['OPLS4', 'OPLS5'], default='OPLS4', help='Specify the forcefield to use. Default: OPLS4.')
    parser.add_argument('--oplsdir', type=str, default=None)
    parser.add_argument('--ffbuilder', action='store_true')
    parser.add_argument('--ffhost', type=str, default='compute-16core-64gb-ondemand')
    parser.add_argument('--md-sim-time', type=float, default=2000)
    parser.add_argument('--fep-sim-time', type=float, default=10000)
    parser.add_argument('--salt', type=float, default=None)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--lig-rest', action='store_true')
    parser.add_argument('--seed', type=int, default=2007)
    args = parser.parse_args()

    # Remove FFB host if no ff is generated
    if not args.ffbuilder:
        ffhost = None
    # File paths from input 
    mae_path = Path(args.mae_file)
    job_path = Path(args.job_name)
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
    write_abfep_job_script(
        args.job_name, 
        args.fep_sim_time, args.md_sim_time, 
        args.host, args.subhost, 
        args.project, args.maxjob, args.retries, 
        salt=args.salt,
        ff=args.ff,
        opls=opls, 
        ffhost=args.ffhost,
        seed=args.seed
    )        
    # Leave the job directory
    os.chdir('..')

