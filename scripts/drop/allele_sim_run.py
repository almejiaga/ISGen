#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 13:08:27 2015

@author: dominic
"""

from __future__ import print_function
import sys,os
sys.path.append(os.path.abspath('../utils/'))
import allele_sim as sim
import build_panel as build
#import matplotlib.pyplot as plt

eprint = lambda *args, **kwargs: print(*args, file=sys.stderr, **kwargs)

## Old plotting
#plt.xlabel('Panel Size')
#plt.ylabel('Coverage')
#plt.plot(x, greedycoverage, label = 'greedy')
#plt.plot(x, randomcoverage, label = 'random')
#plt.xlim(0, x[-1] + 1)
#plt.legend(loc = 4)
#plt.show()

class CLIError(Exception):
    pass

try:
    # simulation options
    iteration_count = None
    input_path = None
    new_file_decision = "new" # Either --new-file, --append, or --overwrite

    # building options
    panel_output_path = "."
    panel_size_list = None
    random_panel_count = 0

    # shared
    allele_simulation_path = "."

    # process type
    simulate = None
    do_build = None

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        n = lambda: sys.argv[i+1]

        if arg == '--simulate':
            simulate = True
        elif arg == '--build':
            do_build = True
        elif arg == '--iterations':
            iteration_count = int(n())
            i += 1
        elif arg == '--input':
            input_path = n()
            i += 1
        elif arg == '--proband-file':
            proband_file = n()
            i += 1
        elif arg == '--sim-output-dir':
            allele_simulation_path = n()
            i += 1
        elif arg == '--new-file':
            new_file_decision = "new"
        elif arg == '--append':
            new_file_decision = "append"
        elif arg == '--overwrite':
            new_file_decision = "overwrite"
        elif arg == '--allele-simulation':
            allele_simulation_path = n()
            i += 1
        elif arg == '--panel-output-dir':
            panel_output_path = n()
            i += 1
        elif arg == '--panel-size-list':
            panel_size_list = map(int, n().split(','))
            i += 1
        elif arg == '--random-points':
            random_panel_count = int(n())
            i += 1
        else:
            raise CLIError("unrecognized command-line option %s" % arg)
        i += 1
except CLIError as e:
    eprint("Command-line interface error:", e)
    sys.exit(1)
except IndexError:
    eprint("Unexpected end of command-line arguments.")
    sys.exit(1)

def check_arg(arg_name, arg):
    if arg is None:
        raise CLIError('missing mandatory argument %s' % arg_name)

if do_build is None and simulate is None:
    eprint("What are you doing with your life?")
    sys.exit(42)

try:
    # check arguments:
    if simulate:
        check_arg('iterations', iteration_count)
        check_arg('input', input_path)

    if (not simulate) and build:
        check_arg('sim-output-dir', allele_simulation_path)

    if do_build:
        check_arg('panel-output-dir', panel_output_path)
        check_arg('panel-size-list', panel_size_list)
        check_arg('random-points', random_panel_count)
		#Must be specified as None if you want all probands to be in
        check_arg('proband-file', proband_file)


    # Do the work.
    if simulate:
        results = sim.run(
                os.path.expanduser(input_path),
                os.path.expanduser(allele_simulation_path),
                iteration_count,
                new_file_decision)
        sim_results_path = results[0]
    else:
        sim_results_path = allele_simulation_path

    if do_build:
        print("Finding allele matches among individuals...")
        		
        panel = build.panel(sim_results_path,proband_file)
        eprint("Done")

        print("Read", panel.Iterations, "simulations from file.")
        greedycoverage, randomcoverage, runtimes = panel.run(panel_output_path,
                                                             panel_size_list,
                                                             random_panel_count)
except CLIError as e:
    eprint(e)
    sys.exit(1)
