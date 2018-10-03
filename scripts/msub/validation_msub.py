import sys, os
sys.path.append(os.path.abspath('control/'))
sys.path.append(os.path.abspath('msub/'))
sys.path.append(os.path.abspath('../scripts'))
import argparse
import numpy as np
import subprocess
import copy
import time, datetime
import msub_tools as mt
import pandas as pd
import run


def load_panel_params(panel_file):
    """
    Loads all panels, along with generating ancestor and regional
    allele frequencies.
    """
    panel_data = pd.DataFrame.from_csv(panel_file, sep='\t')

    return panel_data


def make_timestamped_dir(path):
    """ Creates a timestamped directory in the given path """
    date = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    newdir = os.path.join(os.path.expanduser(path), date)
    mt.make_if_not_exists(newdir)

    return newdir


def write_job_panel(out_dir, anc, panel_inds, region_freqs):
    """ Write panel and real allele frequencies to be validated by the job """
    panel_file = os.path.join(out_dir, str(anc) + '_panel.txt')
    region_freqs_file = os.path.join(out_dir, 'expected_freqs.txt')

    with open(panel_file, 'w') as f:
        for ind in panel_inds:
            f.write(ind + '\n')

    with open(region_freqs_file, 'w') as f:
        f.write('Region\tNum_Alleles/Num_Probands\n')
        for key, value in region_freqs.iteritems():
            f.write(key + ': ' + value + '\n')

    return panel_file, region_freqs_file


def get_freq_params(region_freqs, sample_size=1000):
    """
    Reads allele frequency parameters from panel file, and returns parameters
    in tuple (pop_size, sample_size, sample_count), scaled to provided
    sample size
    """
    global_params = region_freqs['All Probands'].split('/')
    global_count, num_probands = map(int, global_params)

    ## NOTE: Add option for choosing sample size +p3 id:183
    ## Reduce sample size for small pedigrees
    sample_size = np.min([sample_size, num_probands])
    
    global_freq = 1. * global_count / num_probands
    sample_count = int(round(global_freq * sample_size))

    return num_probands, sample_size, sample_count


def main(args):
    ## Load panels from file
    panels = load_panel_params(args.panel_file)

    ## Create argparse namespace with args common to all validation jobs,
    ## and filter out args used to setting up batch execution
    batch_args = ['outpath', 'msub_validation', 'panel_file',
                          'resource_header', 'test', 'hom_het']
    job_args = mt.filter_args(args, batch_args, exclude=True)

    for anc, panel_data in panels.iterrows():
        ## Separate clause from allele frequency data
        region_freqs = panel_data.to_dict()
        panel_inds = region_freqs.pop('Clause').split(',')

        ## Create timestamped directory
        ## TODO: Add a flag where all output files are given exactly +p3 id:167
        ## since this will facilitate rerunning partial jobs
        out_dir = make_timestamped_dir(args.outpath)

        ## Write panel inds and regional frequencies to file
        panel_file, region_freqs_file = write_job_panel(out_dir,
                                                    anc,
                                                    panel_inds,
                                                    region_freqs)

        ## Set hom/het panel files
        if args.hom_het == 'het':
            job_args.hetfile = panel_file
        elif args.hom_het == 'hom':
            job_args.homfile = panel_file

        ## DONE: Read allele frequency from file and pass as arg +p2 id:169
        ## Set observed allele frequency parameters
        freq_params = get_freq_params(region_freqs)
        job_args.freq_params = ','.join([str(x) for x in freq_params])

        ## Set output files if they do not exist
        ## TODO: Use provided filenames (only path used here) +p3 id:168
        if not args.climb_outfile or not os.path.isfile(args.climb_outfile):
                job_args.climb_outfile = os.path.join(out_dir, 'climb_out.h5')

        if not args.control_outfile or not os.path.isfile(args.control_outfile):
                job_args.control_outfile = os.path.join(out_dir,
                                                        'control_out.h5')

        if args.calc_regional_freqs and (not args.regional_freq_outfile or
                            not os.path.isfile(args.regional_freq_outfile)):
                job_args.regional_freq_outfile = os.path.join(out_dir,
                                                        'regional_freqs.h5')

        ## Make directory to store job scripts
        script_outpath = os.path.join(out_dir, 'job_script')
        mt.make_if_not_exists(script_outpath)

        ## Create new filename for script, and write header
        script_name = os.path.join(script_outpath, "script.sh")
        subprocess.call(['cp ' + args.resource_header + ' ' + script_name],
                                                            shell=True)

        ## Command to change diretory to job submission directory
        chdir = ["cd $PBS_O_WORKDIR"]

        ## Convert argparse args into CLI compatible command, and add to kwargs
        command = 'python run.py'
        commands = chdir + [mt.argparse_to_CLI(command, job_args)]
        kwargs = {'outpath': script_outpath}
        kwargs['commands'] = commands

        ## Write commands to resource header file
        with open(os.path.expanduser(script_name), 'a') as f:
            for line in mt.build_job_script(kwargs, append=True):
                f.write(line)

        ## Submit job and wait, so that timestamped directories will be unique
        if not args.test:
            subprocess.call(['qsub ' + script_name], shell=True)
        time.sleep(1)
