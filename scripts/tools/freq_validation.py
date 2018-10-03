import sys, os
module_path = os.path.expanduser('~/project/anc_finder/scripts/')
sys.path.append(module_path)
sys.path.append(module_path + 'climb/')
sys.path.append(module_path + 'control/')
sys.path.append(module_path + 'msub/')
import argparse
import run
import tempfile
import tables
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import numpy as np
import msub_tools as mt


def iter_panels(panels_file):
    with open(panels_file, 'r') as f:
        header = next(f).split('\t')
        regions = [x.strip() for x in header[2:]]
        pop_region_ix = 0
        for x in header:
            if x == 'All Probands':
                break
            pop_region_ix += 1
        params = {}
        for line in f:
            line = line.strip().split('\t')
            anc = int(line[0])
            panel = [int(ind) for ind in line[1].split(',')]
            nobs = np.zeros(len(regions), dtype=int)
            sizes = np.zeros(len(regions), dtype=int)
            for i, r in enumerate(regions):
                ## Offset by 2 to skip anc and panel entries
                nobs[i] = int(line[i+2].split('/')[0])
                sizes[i] = int(line[i+2].split('/')[1])

            freqs = pd.DataFrame()
            freqs['nobs'] = nobs
            freqs['N'] = sizes
            freqs['region'] = regions
            
            yield anc, panel, freqs


def write_panel(panel_inds, outfile=None):
    if outfile is None:
        panel_file = tempfile.NamedTemporaryFile(delete=False)
        for ind in panel_inds:
            panel_file.write(str(ind) + '\n')

        return panel_file.name
    else:
        with open(outfile, 'w') as f:
            for ind in panel_inds:
                f.write(str(ind) + '\n')

        return outfile



def cleanup(fname):
    if os.path.exists(fname):
        os.remove(fname)


def read_freqs(freq_outfile):
    with tables.open_file(freq_outfile) as f:
        freqs = f.root.regional_expected_freqs[:]

    return freqs


def simulate(run_args, qsub_args=None):
    freq_estimate = None
    try:
        if qsub_args is None:
            run.main(run_args)
            freq_estimate = pd.DataFrame(read_freqs(run_args.regional_freq_outfile),
                    columns = ['Region', 'Estimate', 'N'])
        else:
            import subprocess, time
            script_dir = os.path.expanduser('~/project/anc_finder/scripts/')

            command = mt.argparse_to_CLI('python run.py', run_args)
            subprocess.call(['python ' +\
                    script_dir + 'msub/quick_submit.py ' +\
                    qsub_args.header_file + ' ' +\
                    qsub_args.outdir + ' ' +\
                    script_dir + ' ' +\
                    command], shell=True)
            time.sleep(2)
    except:
        cleanup(run_args.climb_outfile)
        cleanup(run_args.control_outfile)
        cleanup(run_args.regional_freq_outfile)
        raise

    return freq_estimate


def main(args):

    hetfile = None
    homfile = None

    panels = iter_panels(os.path.expanduser(args.panelsfile))
    true_counts = []
    estimated_counts = []


    for i, params in enumerate(panels):
        if args.num_panels > 0 and i >= args.num_panels:
            break

        timestamped_dir = mt.make_timestamped_dir(args.outdir)
        climbfile = timestamped_dir + '/climb.h5'
        controlfile = timestamped_dir + '/control.h5'
        freq_outfile = timestamped_dir + '/freq.h5'
        anc, panel, true_freqs = params

        print "Estimating allele frequency for panel", i + 1, \
                "consisting of", panel
        print "True simulated allele frequencies"
        print true_freqs
        true_freqs.to_csv(timestamped_dir + '/true_freqs.txt')
        if args.homhet == 'het':
            hetfile = write_panel(panel, timestamped_dir + '/hets.txt')
        elif args.homhet == 'hom':
            homfile = write_panel(panel, timestamped_dir + '/homs.txt')

        N = true_freqs[true_freqs['region'] == 'All Probands']['N'].values[0]
        nobs = true_freqs[true_freqs['region'] == 'All Probands']['nobs'].values[0]
        sample_size = min(1000, N)
        freq_params = str(N) + ',' + str(sample_size) + ',' + \
                str(int(np.round(1. * sample_size * nobs / N)))
        print "Frequency parameters:", freq_params

        run_args = argparse.Namespace(
                climb_iterations=args.iterations,
                contrib_outfile=os.path.expanduser(args.contribfile),
                pedfile=os.path.expanduser(args.pedfile),
                climb_outfile=climbfile,
                control_outfile=controlfile,
                freq_params=freq_params,
                regional_freq_outfile=freq_outfile,
                hetfile=hetfile,
                homfile=homfile,
                verbose=True,
                regions_to_calculate=None, # Estimates all regions by default
                sim_homs=0.5
                )

        qsub_args = None
        if args.qsub is True:
            qsub_args = argparse.Namespace(
                    header_file=args.resource_header,
                    outdir=timestamped_dir)

        freq_estimate = simulate(run_args, qsub_args)
        if freq_estimate is not None:
            freq_estimate.to_csv(timestamped_dir + '/estimated_freqs.txt')
        #
        # true_counts.append(true_freqs.sort_values('region').values)
        # if freq_estimate is not None:
        #     estimated_counts.append(freq_estimate.sort_values('Region').values)

        # cleanup(climbfile)
        # cleanup(controlfile)
        # cleanup(freq_outfile)
        #
    # all_true_counts = np.vstack(true_counts)[:, 0]
    # all_est_counts = np.vstack(estimated_counts)[:, 1].astype(float)
    # print "True freqs"
    # print all_true_counts
    # print "Estimated freqs"
    # print all_est_counts
    #
    # outdir = os.path.expanduser(args.outdir)
    # with open(outdir + '/freq_results.txt', 'w') as f:
    #     f.write('True\tEstimated\n')
    #     for true, est in zip(all_true_counts, all_est_counts):
    #         f.write(str(true) + '\t' + str(est) + '\n')
    #
    # plt.scatter(all_true_counts, all_est_counts)
    # plt.title('True vs Estimated Global Allele Frequencies')
    # plt.xlabel('True Allele Frequency')
    # plt.ylabel('Estimated Allele Frequency')
    # plt.savefig(outdir + '/freq_results.png')
    #
    # return all_true_counts, all_est_counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pedfile", metavar="|",
            help="Pedigree for which to validate allele frequency estimates",
            required=True)
    parser.add_argument("-c", "--contribfile", metavar="|",
            help="File containing results of pedigree-wide allele dropping" +\
                    "simulations", required=True)
    parser.add_argument("-a", "--panelsfile", metavar="|",
            help="File containing simulated patient panels, along with " +\
                    "associated regional allele frequencies", required=True)
    parser.add_argument("-i", "--iterations", metavar="|",
            help="Number of climbing iterations to perform", type=int,
            required=True)
    parser.add_argument("-n", "--num_panels", metavar="|",
            help="Number of panels for which to estimate allele frequencies",
            type=int, default=-1)
    parser.add_argument("-o", "--outdir", metavar="|",
            help="Directory in which to store output", default='')
    parser.add_argument("-H", "--homhet", metavar="|",
            help="Format ['hom' (default) | 'het'] - specify whether " +\
                    "provided panels are homozygotes or heterozygotes",
            default='hom')
    parser.add_argument('-q', '--qsub',
            help='Flag to submit each panel as a separate job to qsub',
            action='store_true')
    parser.add_argument('-r', '--resource_header', metavar='|',
            help='Header containing resource allocations for qsub jobs')
    args = parser.parse_args()

    if args.qsub is True:
        assert args.resource_header is not None, "Must provide resource header"

    main(args)


