from __future__ import division
import tables
import sys, os
sys.path.append('../bootstrap')
import bootstrap
import traceback
import numpy as np
import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from scipy.stats import beta 
from collections import Counter


def get_posterior(sample_df):
    result = bootstrap.Bootstrap.alpha_posterior(sample_df).transpose()
    result.reset_index(level=0)

    return pd.DataFrame(result)


def label_gen_anc(sym_means, sym_gen_anc):
    """
    Adds a column with a binary indicator whether the ancestor generated
    the panel or not
    """
    results = pd.DataFrame(sym_means['posterior'])
    results['generator'] = np.zeros(results.shape[0])

    try:
        results.ix[sym_gen_anc]['generator'] = 1.
    except KeyError:
        ## Add the generating ancestor if they are not present in the results
        gen_anc_df = pd.DataFrame.from_dict({'posterior':{sym_gen_anc:0},
                                         'generator':{sym_gen_anc:1}})
        results = results.append(gen_anc_df)

    return results


def get_mean_accuracy(all_means, nbins=10):
    """
    Bins ancestors according to mean bootstrapped posterior probability,
    and then returns the mean accuracy for each bin
    """
    ## Add a columns of bin assignments
    # bins = np.linspace(0, all_means['posterior'].max(), nbins)
    bins = np.linspace(0, 1, nbins)
    all_means['bin'] = np.digitize(all_means['posterior'], bins)

    ## Add upper bound to right-most bin
    all_means.replace(to_replace={'bin':{nbins: nbins-1}}, inplace=True)

    ## Bin ancestors by mean bootstrapped probability, adding columns for
    ## whether they were the true generating ancestor, and the number of
    ## ancestors in each bin
    bin_count = lambda x: len(x)
    binned = all_means[['generator', 'bin']].pivot_table(index='bin',
                            aggfunc=[np.mean, bin_count], fill_value=0)
    binned.columns = [['observed_prob', 'bin_count']]
    binned['n_successes'] = binned['observed_prob'].values * \
            binned['bin_count'].values

    ## Estimate means and confidence intervals as sampling from a binomial
    ## distribution, with a uniform prior on success rates - Done using
    ## a beta distribution
    binned['alpha'] = binned['n_successes'] + 1
    binned['beta'] = binned['bin_count'].values - binned['n_successes'].values + 1
    beta_mean = lambda row: beta.mean(float(row['alpha']), float(row['beta']))
    binned['posterior_mean'] = binned.apply(beta_mean, axis=1)

    ## Add confidence intercals
    beta_025CI = lambda row: beta.ppf(0.025, float(row['alpha']), float(row['beta']))
    beta_975CI = lambda row: beta.ppf(0.975, float(row['alpha']), float(row['beta']))
    binned['CI2.5'] = binned.apply(beta_025CI, axis=1)
    binned['CI97.5'] = binned.apply(beta_975CI, axis=1)

    ## Convert to values relative to mean, to fit plotting convention
    binned['CI2.5'] = binned['posterior_mean'].values - binned['CI2.5'].values
    binned['CI97.5'] = binned['CI97.5'].values - binned['posterior_mean'].values

    ## Add column with bin centre for plotting
    binned['bin_centre'] = all_means[['posterior', 'bin']].groupby('bin').mean()

    return binned


def plot_accuracy(binned, outfiles):
    colours = itertools.cycle(sns.color_palette())
    point_colour = next(colours)
    line_colour = next(colours)
    ## Plot each point individually so we can have variable point sizes
    x_vals = binned['bin_centre'].values
    y_vals = binned['posterior_mean'].values
    low_err = binned['CI2.5'].values
    high_err = binned['CI97.5'].values
    bin_counts = binned['bin_count'].values
    lim = np.max([np.max(x_vals), np.max(y_vals + high_err)])

    matplotlib.rcParams.update({'font.size': 16})

    for x, y, l, h, c in zip(x_vals, y_vals, low_err, high_err, bin_counts):
        plt.errorbar(x, y, yerr=([l], [h]), capsize=4, fmt='o',
                color=point_colour, markersize=np.log(c) / np.log(1.4))
    # plt.title('Validation of inferred posterior probability')
    plt.ylabel('Proportion of True Ancestors')
    plt.xlabel('Inferred Posterior Probability of Ancestor')
    margin=0.05
    plt.xlim((0-margin, lim+margin))
    plt.ylim((0-margin, lim+margin))
    plt.tight_layout()
    plt.plot([0, lim+margin], [0, lim+margin], color=line_colour)
    for outfile in outfiles:
        plt.savefig(outfile)
    plt.clf()
    import IPython; IPython.embed()


def read_all_runs(hdf5_file):
    """
    Loads all runs from nodes in the provided hdf5 file.
    """
    with pd.HDFStore(hdf5_file) as store:

        for key in store.keys():
            df = store.get(key)

            ## Load ancestor (up to symmetry) who generated the panel for
            ## this validaiton sim
            sym_gen_anc = store.get_storer(key).attrs.sym_gen_anc

            if df.shape[0] == 300000:
                yield df, sym_gen_anc
            else:
                print "Wrong shape!", df.shape


def integral_accuracy(binned_probs):
    """
    Calculates the area under the curve delimited by the estimated
    ancestor posterior probabilities. A value of 0.5 is consistent
    with perfectly calibrated estimates.
    """
    ##TODO: Needs to be updated - bin_width is no longet constant
    x_vals = binned_probs['bin_centre'].values
    y_vals = binned_probs['posterior', 'generator'].values

    bin_width = x_vals[1] - x_vals[0]
    integral = np.sum(np.multiply(bin_width, y_vals))

    return integral


def main(args):
    if args.n_samples is not None:
        sample_sizes = map(int, args.n_samples.split(','))
    else:
        df, _ = next(read_all_runs(args.liks_file))
        sample_sizes = [df.shape[0]]

    all_accuracy = np.zeros((len(sample_sizes), args.iterations))
    
    for i, sample_size in enumerate(sample_sizes):
        for j in range(args.iterations):
            ## Load data from CLI args
            print "Loading likelihoods..."
            data = read_all_runs(args.liks_file)

            ## Calculate the likelihood of each ancestor
            all_liks = []
            for df, sym_gen_anc in data:

                ## Subsample data, or skip if not enough iterations
                ## are present
                try:
                    df = df.sample(sample_size, replace=False)
                except ValueError:
                    continue

                ## Append indicator column denoting if ancestor generated
                ## panel
                anc_liks = get_posterior(df)
                anc_liks.columns = ['posterior']
                all_liks.append(label_gen_anc(anc_liks, sym_gen_anc))

            ## Concatenate all results and write to file
            all_results = pd.concat(all_liks)
            binned = get_mean_accuracy(all_results, nbins=args.num_bins)
            path, ext = os.path.splitext(args.outfile)
            label = str(sample_size) + '_' + str(j)

            ## Write raw accuracy results
            raw_outfile = path + '_raw_' + label + ext
            binned.to_csv(raw_outfile)

            ## Plot inferred vs empirical probabilities
            plot_outfiles = []
            plot_outfiles.append(path + '_plot_' + label + '.png')
            plot_outfiles.append(path + '_plot_' + label + '.pdf')
            plot_outfiles.append(path + '_plot_' + label + '.svg')
            plot_outfiles.append(path + '_plot_' + label + '.ps')
            plot_accuracy(binned, plot_outfiles)

            ## Evaluate integral accuracy
            # accuracy = integral_accuracy(binned)
            # all_accuracy[i, j] = accuracy
            # print "Accuracy for", sample_size, "is", accuracy

    accuracy_df = pd.DataFrame(all_accuracy.transpose())
    accuracy_df.columns = sample_sizes
    print "Average accuracy:", accuracy_df.mean(axis=0)
    accuracy_df.to_csv(args.outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--num-bins", metavar='|',
                        help="Number of evenly-spaced bins to use when " +\
                            "partitioning posterior probabilities",
                        type=int, default=11)
    parser.add_argument("-s", "--n_samples", metavar='|',
                        help="Set of number of rows to sample from each " +\
                        "simulation, in format n1,n2,n3,...,n_last",
                        default=None)
    parser.add_argument("-i", "--iterations", metavar='|',
                        help="Number of times to calculate accuracy of " +\
                        " specified number of samples",
                        default=1, type=int)

    ## DONE: Some of these are required only if file does not exist +p2 id:119
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-l", "--liks_file", metavar='|',
                        help="File containing all ancs and corresponding " +\
                        "likelihoods for a batch of validation simulations",
                        required=True)
    requiredNamed.add_argument("-o", "--outfile", metavar='|',
                        help="File to output accuracy for each likelihood" +\
                        " bin, and plot of observed vs expected accuracy",
                        required=True)

    args = parser.parse_args()

    ## Args for paper:
    # -l ~/project/anc_finder/results/validation_pedEx/all_liks.h5
    # -n 7

    ## This is set here so that we can check whether we are calling from
    ## the CLI or another script above
    args.from_file = True
    args.CIs = None

    main(args)
