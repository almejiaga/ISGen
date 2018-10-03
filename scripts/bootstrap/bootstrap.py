from __future__ import division
import numpy as np
import pandas as pd
import sys, os
import traceback
from collections import defaultdict, Counter
import argparse
import tables
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn


class Bootstrap:
    def __init__(self, climb_file, control_file, symfile=None,
                 jackknife_subsamples=None):
        self.jackknife_subsamples = jackknife_subsamples

        ## Load ancs and total (climb * control) likelihoods
        with tables.open_file(climb_file, 'r') as f:
            self.ancs = f.root.ancs[:]
        with tables.open_file(control_file, 'r') as f:
            if 'tot_hap_liks' in [node.name for node in f.list_nodes(f.root)]:
                liks = f.root.tot_hap_liks[:]
                print "Using haplotype likelihoods"
            else:
                liks = f.root.tot_liks[:]
                print "No haplotype likelihoods found"

            self.liks = np.asarray(map(lambda x: np.exp2(float(x)), liks))

        self.likdata = pd.DataFrame(zip(self.ancs, self.liks))
        self.likdata.columns = ["Ancs", "Liks"]

        if symfile is not None:
            symdat = pd.DataFrame.from_csv(symfile)
            print len(set(self.ancs)), "original ancestors"

            ## Split into original ancs and their symmetrical replacements
            symdat.columns = ['Symmetrical']
            symdat['Original'] = symdat.index
            symdat.reset_index(level=0, inplace=True)

            ## Collapse symmetrical ancestors
            self.likdata.replace(to_replace=symdat['Original'].tolist(),
                                 value=symdat['Symmetrical'].tolist(),
                                 inplace=True)

            print len(set(self.likdata['Ancs'])), "after collapsing symmetries"
            self.symdat = symdat
        else:
            self.symdat = None


    def bootstrap(self, iterations, method=None, verbose=True):
        ## Bootstraps output of 'self.likfile', transforming it using 'method'
        if method is None:
            method = Bootstrap.alpha_posterior

        ## Take jackknife subsample if specified
        if self.jackknife_subsamples is None:
            boot_data = self.likdata
        else:
            boot_data = self.likdata.sample(n=self.jackknife_subsamples,
                                       replace=False)

        ## Store results in a list to concatenate later
        bootresults = []

        ## Perform bootstrap
        for i in xrange(iterations):
            if verbose:
                if i % 100 == 0: print "Iteration", i, "of", iterations
            resamp = boot_data.sample(frac=1, replace=True)
            result = method(resamp).transpose()
            result.reset_index(level=0)
            bootresults.append(result)

        bootresults = pd.DataFrame.from_dict(map(dict,bootresults))

        return bootresults


    @staticmethod
    def confidence_intervals(bootresults):
        intervals = bootresults.describe([0.025,0.975]).transpose()
        CIs = intervals[['2.5%', 'mean', '97.5%']]

        ## If an ancestor is not present in a resample, they get probability 0
        CIs = CIs.fillna(0)

        return CIs


    @staticmethod
    def save_boot(bootresults, outfile, format='confidence', zerofill=True):
        if zerofill is True:
            bootresults = bootresults.fillna(0)

        if format == 'raw':
            bootresults.to_csv(outfile)
        elif format == 'confidence':
            CIs = Bootstrap.confidence_intervals(bootresults)
            path, ext = os.path.splitext(outfile)
            plot_outfile = path + '_plot.png'
            CIs.to_csv(outfile)
            Bootstrap.plot_CIs(CIs, plot_outfile)
        else:
            print "Unrecognized format:", format
            print "Aborting..."


    @staticmethod
    def sort_liks(sample_pd):
        group = sample_pd.groupby(['Ancs'], sort=False)
        sort_by_anc = group.sum() / (1. * len(sample_pd))

        return sort_by_anc


    @staticmethod
    def alpha_posterior(sample_pd):
        alphas = Bootstrap.sort_liks(sample_pd)
        tot_lik = alphas.sum()['Liks']

        posterior = alphas['Liks'] / tot_lik

        return posterior


    @staticmethod
    def plot_CIs(CIs, outfile):
        fig, ax = plt.subplots()

        ## Sort values
        CIs = CIs.sort_values(by='mean', ascending=False)

        ## Get means and convert confidence intervals into relative errors
        ## so they can be plotted using plt.errorbar
        means = CIs['mean'].values
        x = np.arange(len(means))
        low_err = means - CIs['2.5%'].values
        high_err = CIs['97.5%'].values - means

        ## Plot and save
        ax.errorbar(x, means, yerr=[low_err, high_err], capthick=2, fmt='o')
        plt.xticks(x, CIs.index, rotation='vertical')
        ax.set_title('Bootstrapped posterior probabilities with 95% \n' +\
                     'confidence intervals', fontsize=18)
        ax.set_ylabel('Posterior Probability', fontsize=18)
        ax.set_xlabel('Ancestor', fontsize=18)
        ax.margins(0.02)
        fig.set_tight_layout(True)
        plt.savefig(outfile)


def main(args):
    ## Set arguments common to both 'batch' and 'single' operation
    iterations = args.iterations
    verbose = not args.quiet

    ## If processing single file, put arguments into lists with one item,
    ## so we can still iterate through them
    if args.subparser_name == 'single':
        climb_files = [os.path.expanduser(args.climb_file)]
        control_files = [os.path.expanduser(args.control_file)]

        if args.symmetry_file is not None:
            symmetry_files = [os.path.expanduser(args.symmetry_file)]
        else:
            symmetry_files = [None]

        if args.raw_outfile is not None:
            raw_outfiles = [os.path.expanduser(args.raw_outfile)]
        else:
            raw_outfiles = [None]

        if args.conf_outfile is not None:
            conf_outfiles = [os.path.expanduser(args.conf_outfile)]
        else:
            conf_outfiles = [None]

    ## Format batch options so we can iterate through them
    if args.subparser_name == 'batch':
        climb_files = map(os.path.expanduser,
                    [line.strip() for line in open(args.climb_files, 'rU')])
        control_files = map(os.path.expanduser,
                    [line.strip() for line in open(args.control_files, 'rU')])

        if args.symmetry_file is not None:
            symmetry_files = map(os.path.expanduser,
                    [line.strip() for line in open(args.symmetry_file, 'rU')])
        else:
            symmetry_files = [None for i in range(len(climb_files))]

        if args.raw_outfile is not None:
            raw_outfiles = map(os.path.expanduser,
                        [line.strip() for line in open(args.raw_outfile, 'rU')])
        else:
            raw_outfiles = [None for i in range(len(climb_files))]

        if args.conf_outfile is not None:
            conf_outfiles = map(os.path.expanduser,
                    [line.strip() for line in open(args.conf_outfile, 'rU')])
        else:
            conf_outfiles = [None for i in range(len(climb_files))]

    ## Perform the bootstrapping, storing results so they can be returned
    ## to another script if needed
    all_CIs = []
    for climb_file, control_file, symfile, raw_outfile, conf_outfile in \
                zip(climb_files, control_files, symmetry_files,
                    raw_outfiles, conf_outfiles):
        try:
            print "Loading likelihoods..."
            B = Bootstrap(climb_file, control_file, symfile=symfile,
                          jackknife_subsamples=args.jackknife_subsamples)
            print "Done.\n"

            print "Bootstrapping..."
            bootresults = B.bootstrap(iterations, verbose=verbose)
            all_CIs.append(B.confidence_intervals(bootresults))
            print "Done.\n"

            if raw_outfile is not None:
                print "Saving raw output to", raw_outfile
                B.save_boot(bootresults, outfile=raw_outfile, format='raw')
                print "Done.\n"

            if conf_outfile is not None:
                print "Saving confidence intervals to", conf_outfile
                B.save_boot(bootresults, outfile=conf_outfile,
                            format='confidence')
                print "Done.\n"

            if args.raw_outfile is None and args.conf_outfile is None \
                                                               and verbose:
                print "Bootstrap Confidence Intervals:"
                print B.confidence_intervals(bootresults)

            print "Done.\n\n"

        except KeyboardInterrupt:
            raise
        except:
            print "Error!"
            with open('err.log', 'a') as f:
                f.write('-'*60)
                f.write('Error processing the following files:\n')
                f.write(str(climb_file) + '\n')
                f.write(str(control_file) + '\n')
                f.write(str(symfile) + '\n')
                f.write(str(raw_outfile) + '\n')
                f.write(str(conf_outfile) + '\n')
                traceback.print_exc(file=f)
            continue

    return all_CIs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    sp = parser.add_subparsers(dest='subparser_name')

    ## Arguments for bootstrapping single likelihood file
    single_parser = sp.add_parser("single", help="Bootstrap single file")
    single_req = single_parser.add_argument_group("Required named arguments")
    single_req.add_argument("-l", "--climb-file", metavar='',
                        help="File containing the results of allele " +\
                             "climbing simulations",
                        required=True)
    single_req.add_argument("-o", "--control-file", metavar='',
                        help="File containing control likelihoods for " +\
                             "the trees simulated in --climb-file",
                        required=True)
    single_req.add_argument("-i", "--iterations", metavar='',
                        help="Number of bootstrap iterations to perform",
                        required=True, type=int)

    single_parser.add_argument("-s", "--symmetry-file", metavar='',
                        help="File containing individuals who are expected " +\
                        "to have symmetrical likelihoods, based on a " +\
                        "set of probands")
    single_parser.add_argument("-r", "--raw-outfile", metavar='',
                        help="File to output values of all bootstrap resamples")
    single_parser.add_argument("-c", "--conf-outfile", metavar='',
                        help="File to output mean value of bootstrap " +\
                        "resamples, along with 95 percent confidence intervals")
    single_parser.add_argument("-q", "--quiet",
                        help="Suppress printing number of iterations completed",
                        action='store_true')
    single_parser.add_argument("-j", "--jackknife_subsamples", metavar='',
                        help="Number of jackknife subsamples to extract " +\
                        "from the data", type=int)

    ## Arguments for bootstrapping batch of likelihood file
    batch_parser = sp.add_parser("batch", help="Bootstrap batch of files")
    batch_req = batch_parser.add_argument_group("Required named arguments")
    batch_req.add_argument("-l", "--climb-files", metavar='',
                        help="File listing output files of allele " +\
                             "climbing simulations",
                        required=True)
    batch_req.add_argument("-o", "--control-files", metavar='',
                        help="File listing control likelihoods for " +\
                             "the trees simulated in --climb-files",
                        required=True)
    batch_req.add_argument("-i", "--iterations", metavar='',
                        help="Number of bootstrap iterations to perform",
                        required=True, type=int)

    batch_parser.add_argument("-s", "--symmetry-file", metavar='',
                    help="File listing symmetry files to use with " +\
                    "bootstrapping")
    batch_parser.add_argument("-r", "--raw-outfile", metavar='',
                        help="File containing paths to output values of all" +\
                        " bootstrap resamples")
    batch_parser.add_argument("-c", "--conf-outfile", metavar='',
                        help="File containing paths to output mean value of" +\
                        " bootstrap resamples, along with 95 percent" +\
                        " confidence intervals")
    batch_parser.add_argument("-q", "--quiet",
                        help="Suppress printing number of iterations completed",
                        action='store_true')
    batch_parser.add_argument("-j", "--jackknife_subsamples", metavar='',
                        help="Number of jackknife subsamples to extract " +\
                        "from the data", type=int)

    args = parser.parse_args()

    main(args)
