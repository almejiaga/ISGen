import sys, os
sys.path.append('../control/')
import subprocess
import copy
import time, datetime
import argparse


def make_if_not_exists(path):
    """
    Creates 'path' unless it exists, in which case do nothing.
    """
    try:
        os.makedirs(path)
    except OSError, e:
        ## If error raised because path exists, do nothing
        if e.errno != 17:
            raise


def make_timestamped_dir(path):
    """ Creates a timestamped directory in the given path """
    date = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    newdir = os.path.join(os.path.expanduser(path), date)
    make_if_not_exists(newdir)

    return newdir


def build_job_script(kwargs, append=False):
    """ Writes args to a script which can be submitted via msub """
    ## IDEA: Could also provide a file with header to be used +p3 id:122
    lines = []
    
    if append is False:
        line = "#!/bin/bash\n\n"
        lines.append(line)

    ## First write resource allocations
    if 'project' in kwargs:
        line = "#PBS -A " + kwargs['project'] + "\n"
        lines.append(line)

    if 'walltime' in kwargs:
        line = "#PBS -l walltime=" + kwargs['walltime'] + "\n"
        lines.append(line)

    if 'nodes' in kwargs:
        assert 'ppn' in kwargs, "Must specify processors per node"
        line = "#PBS -l nodes=" + kwargs['nodes'] +\
                ":ppn=" + kwargs['ppn'] + "\n"
        lines.append(line)

    if 'queue' in kwargs:
        line = "#PBS -q " + kwargs['queue'] + "\n"
        lines.append(line)

    if 'name' in kwargs:
        line = "#PBS -N " + kwargs['name'] + "\n"
        lines.append(line)

    if 'outpath' in kwargs:
        line1 = "#PBS -o " + kwargs['outpath'] + "/\n"
        line2 = "#PBS -e " + kwargs['outpath'] + "/\n"
        lines.extend([line1, line2])

    if 'dependencies' in kwargs:
        line = "\n#PBS -W x=depend=afterok"
        for job_name in kwargs['dependencies']:
            line += ":JOBNAME." + job_name
        line += '\n'
        lines.append(line)

    ## Now write commands
    lines.append('\n')
    for command in kwargs['commands']:
        lines.append(command + '\n')

    return lines


def argparse_to_CLI(command, args):
    """ Converts an argparse Namespace to a single-line CLI command """
    for arg, value in vars(args).iteritems():
        ## Don't add args that are flagged as False
        if type(value) is not bool:
            command += " --" + arg + " " + str(value)
        elif value is not False:
            command += " --" + arg

    return command


def filter_args(args, include_args, exclude=False):
    """
    Copies args which are present in 'include_args' and discards the rest.
    If exclude=True, copies args *not* present in include_args.
    """
    args_dict = vars(copy.copy(args))
    filtered_args = {key: value for key, value in args_dict.iteritems()
                            if key in include_args and value is not None
                            and value != "None"}

    new_args = argparse.Namespace()

    ## Add either the filtered args or their complement depending on
    ## exclude flag
    if exclude is False:
        vars(new_args).update(filtered_args)
    else:
        vars(new_args).update({key:args_dict[key] for key, value in
                            args_dict.iteritems() if key not in filtered_args
                            and value is not None and value != "None"})

    print vars(new_args)

    return new_args
