import sys, os
sys.path.append('../control/')
import argparse
import time, datetime
import subprocess
import ind_contrib
import msub_tools as mt


def main(args):
    ## Select resource allocation args to be writted in header of script
    resource_args = ['project', 'walltime', 'nodes', 'ppn']
    resources = vars(mt.filter_args(args, resource_args))

    ## Filter out args to be passed to CLI in job script
    include_args = ['regionfile', 'delimiter', 'verbose',
                    'pedfile']
    job_args = mt.filter_args(args, include_args)

    ## Make directory to store output of each single job
    basepath = os.path.dirname(os.path.expanduser(args.outfile))
    job_outpath = os.path.join(basepath, 'part')
    mt.make_if_not_exists(job_outpath)

    ## Make directory to store job scripts
    script_outpath = os.path.join(basepath, 'job_scripts')
    mt.make_if_not_exists(script_outpath)

    ## Add job-specific args
    njobs = int(args.iterations / args.iter_per_job)
    job_args.iterations = args.iter_per_job

    ## Command to change diretory to job submission directory
    chdir = ["cd $PBS_O_WORKDIR"]

    ## Store job names so we can wait for them to finish before combining
    job_names = []
    job_outfiles = []
    for i in range(njobs):
        ## Create new filename for script
        script_name = os.path.join(script_outpath, "script" + str(i) + ".sh")

        ## Create new filename for job output
        file, ext = os.path.splitext(os.path.basename(args.outfile))
        new_outfile = file + "_" + str(i) + ext + ".part"
        job_args.outfile = os.path.join(job_outpath, new_outfile)
        job_outfiles.append(job_args.outfile)

        ## Create timestamped job name
        date = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        job_name = 'dnelson_job' + str(i) + "-" + date
        job_names.append(job_name)

        ## Convert argparse args into CLI compatible command, and add to kwargs
        command = 'python ../control/ind_contrib.py'
        commands = chdir + [mt.argparse_to_CLI(command, job_args)]
        kwargs = {'name': job_name}
        kwargs.update(resources)
        kwargs['commands'] = commands
        kwargs['outpath'] = script_outpath

        ## Write job script
        with open(os.path.expanduser(script_name), 'w') as f:
            for line in mt.build_job_script(kwargs):
                f.write(line)

        ## Submit job
        if not args.test:
            subprocess.call(['qsub ' + script_name], shell=True)
            time.sleep(2)

    ## Write paths of partial output files
    outpaths_file = os.path.join(job_outpath, "outpaths.txt")
    with open(outpaths_file, 'w') as f:
        for out_file in job_outfiles:
            f.write(out_file + '\n')

    ## Select resource allocation args to be written in header of script
    combine_resource_names = ['project']
    combine_resources = vars(mt.filter_args(args, combine_resource_names))

    ## Allocate resources for combine job
    kwargs = {'name': 'dnelson_combine-' + date,
                'walltime': '10:00:00',
                'nodes': '1',
                'ppn': '5'}
    kwargs.update(combine_resources)

    ## Get outfile from CLI args, and set file to provide paths of partial
    ## output files
    combine_args = mt.filter_args(args, ['outfile'])
    combine_args.outpaths_file = outpaths_file

    ## Write commands for combine script
    command = 'python combine_contribs.py'
    kwargs['commands'] = chdir + [mt.argparse_to_CLI(command, combine_args)]

    ## Combine job depends on the partial simulations having been completed
    kwargs['dependencies'] = job_names
    kwargs['outpath'] = script_outpath

    ## Write script to combine output
    script_name = os.path.join(script_outpath, "combine_script.sh")
    with open(os.path.expanduser(script_name), 'w') as f:
        for line in mt.build_job_script(kwargs):
            f.write(line)

    ## Submit job
    ## TODO: Create job to cancel this job on afternotok +p2 id:124
    if not args.test:
        subprocess.call(['qsub ' + script_name], shell=True)
        time.sleep(2)


if __name__ == '__main__':
    # DONE: Explain that partial files will be stored in current dir +p3 id:120
    ## TODO: Clarify required/optional args id:123
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--regionfile", metavar='|',
                        help="File listing the region of each proband, " +\
                            "with no header")
    parser.add_argument("-m", "--delimiter", metavar='|',
                        help="Delimiter used in --regionfile. If space, " +\
                        "surround by double quotes. Default is tab.")
    parser.add_argument("-v", "--verbose", action="store_true")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-p", "--pedfile", metavar='|',
                        help="Pedigree file",
                        required=True)
    requiredNamed.add_argument("-o", "--outfile", metavar='|',
                        help="File to store simulation output in hdf5 format",
                        required=True)
    requiredNamed.add_argument("-i", "--iterations", metavar='|',
                        help="Number of allele drops to perform from each " +\
                            "individual in the pedigree",
                        type=int, required=True)
    requiredNamed.add_argument("-j", "--iter-per-job", metavar='|',
                        help="Number of iterations per job. Number of jobs" +\
                        " is then determined by iterations / iter-per-job",
                        type=int, required=True)

    jobRequired = parser.add_argument_group('Job submission arguments')
    jobRequired.add_argument("-P", "--project", metavar='|',
                        help="Project to use for resource allocation")
    jobRequired.add_argument("-w", "--walltime", metavar='|',
                        help="Walltime to allocate to each partial job")
    jobRequired.add_argument("-n", "--nodes", metavar='|',
                        help="Number of nodes to request per partial job")
    jobRequired.add_argument("-N", "--ppn", metavar='|',
                        help="Processors requested per node")
    jobRequired.add_argument("-q", "--queue", metavar='|',
                        help="Queue to which to submit jobs")
    jobRequired.add_argument("-t", "--test",
                        help="Write scripts but do not submit them",
                        action='store_true')

    args = parser.parse_args()

    main(args)
