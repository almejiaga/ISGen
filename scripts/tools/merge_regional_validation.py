import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import tables
import pandas as pd
import seaborn
import argparse
import sys, os


def read_freqs(freq_outfile):
    with tables.open_file(freq_outfile) as f:
        freqs = f.root.regional_expected_freqs[:]

    freq_df = pd.DataFrame(freqs, columns=['Region', 'Estimate', 'N'])
    freq_df['Estimate'] = freq_df['Estimate'].astype(float)
    freq_df['N'] = freq_df['N'].astype(float)

    return freq_df


def plot_freq_accuracy(true_freqs, est_freqs, plot_out_file):
    plt.plot([0, 1], [0, 1], '--', color='lightgrey')
    plt.hist2d(true_freqs, est_freqs,
               bins=np.arange(0, 0.01, 0.0005), norm=LogNorm())
    plt.title('True vs Estimated Regional Allele Frequencies')
    plt.xlabel('True Allele Frequency')
    plt.ylabel('Estimated Allele Frequency')
    max_val = np.max([np.max(true_freqs), np.max(est_freqs)])
    plt.xlim(0, 0.01)
    plt.ylim(0, 0.01)
    # plt.xlim(0, max_val)
    # plt.ylim(0, max_val)
    plt.colorbar()
    
    plt.savefig(plot_out_file)
    plt.clf()


def load_data(freq_file, climb_file, control_file):
    df = pd.DataFrame()
    with tables.open_file(freq_file, 'r') as f:
        i = 0
        for node in f.walk_nodes(f.root):
            i += 1
            if i == 1:
                ## Skip the root node - no data
                continue
            if node.name == 'regional_expected_freqs':
                ## Skip combined expectations - we only
                ## want per-tree estimates
                continue

            df[node.name] = node[:]
                
    with tables.open_file(climb_file, 'r') as f:
        df['climb_liks'] = f.root.liks[:]
                
    with tables.open_file(control_file, 'r') as f:
        df['control_liks'] = f.root.control_liks[:]

    return df


def get_expectations(true_freqs, est_df, tree_weight=True):
    """
    Calculates the expected allele frequencies of each region, with two
    options for calibrating results to the observed global frequency.

    tree_weight=True weights each tree by global_freq / E[global_freq | tree]
    tree_weight=False weights all trees by global_freq / E[global_freq]
    """

    if tree_weight is False:
        est_df['weight'] = np.ones(len(est_df['All Probands']))
    else:
        est_df['weight'] = true_freqs['All Probands']['nobs'] \
                                / est_df['All Probands']

    regions = true_freqs.columns.values
    lik_prod = np.exp2(est_df['climb_liks'] + est_df['control_liks'])
    lik_prod_sum = lik_prod.sum()
    true_vals = np.zeros(len(regions))
    est_vals = np.zeros(len(regions))

    global_weight = None
    for j, r in enumerate(regions):
        region_mean = est_df[r].values

        N = true_freqs[r]['N']
        true_vals[j] = 1. * true_freqs[r]['nobs'] / N
        est_vals[j] = (region_mean * est_df['weight'] * lik_prod).sum() \
                        / lik_prod_sum / N
        if r == 'All Probands':
            global_weight = 1. * true_vals[j] / est_vals[j]

    if tree_weight is False:
        assert global_weight is not None
        est_vals = est_vals * global_weight

    return true_vals, est_vals


def main(args):
    outpath = os.path.expanduser(args.outpath)
    true_freqs_files = [line.strip() \
                        for line in open(args.true_freqs_files, 'r')]
    est_freqs_files = [line.strip() \
                        for line in open(args.est_freqs_files, 'r')]
    climb_files = [line.strip() \
                        for line in open(args.climb_files, 'r')]
    control_files = [line.strip() \
                        for line in open(args.control_files, 'r')]
    assert len(true_freqs_files) == len(est_freqs_files)
    plot_out_files = [outpath + '/plot_' + str(i) + '.png' \
                        for i in range(len(est_freqs_files))]

    true_df = pd.DataFrame()
    est_df = pd.DataFrame()
    for i in range(len(true_freqs_files)):
        print "Panel", i
        global_true = None
        global_est = None
        truefile = true_freqs_files[i]
        freqfile = est_freqs_files[i]
        climbfile = climb_files[i]
        controlfile = control_files[i]
        outfile = plot_out_files[i]

        df = load_data(freqfile, climbfile, controlfile)
        true_freqs = pd.DataFrame.from_csv(truefile)
        true_freqs = true_freqs.pivot_table(columns='region')
        regions = true_freqs.columns.values

        true_vals, est_vals = get_expectations(true_freqs, df,
                                               tree_weight=True)
            
        true_df = pd.concat([true_df, pd.DataFrame([true_vals], columns=regions)])
        est_df = pd.concat([est_df, pd.DataFrame([est_vals], columns=regions)])
                
        # inflation = global_est / global_true
        # print("Inflation factor:", inflation)
        ## Per-panel plots can be made here
        # plot_freq_accuracy(true, est, inflation, plot_out_files[0])

    mse = (true_df - est_df) ** 2

    ## Save plot built in plotting function
    all_runs_plot_file = os.path.join(outpath, 'all_runs.png')
    plot_freq_accuracy(true_df.values.ravel(), est_df.values.ravel(),
            all_runs_plot_file)
    mse.to_csv(args.outpath + '/mse.csv')
    true_df.to_csv(args.outpath + '/true_vals.csv')
    est_df.to_csv(args.outpath + '/est_vals.csv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("true_freqs_files", 
                        help="File listing paths to true allele frequencies")
    parser.add_argument("est_freqs_files", 
                        help="File listing paths to estimated allele frequencies")
    parser.add_argument("climb_files",
                        help="File listing paths to climbing likelihood files")
    parser.add_argument("control_files",
                        help="File listing paths to control likelihood files")
    parser.add_argument("outpath", help="Path to output plots")

    args = parser.parse_args()
    main(args)
