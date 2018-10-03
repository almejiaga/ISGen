import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import tables
import pandas as pd
import scipy.sparse
# import seaborn as sns; sns.set()
import argparse
from scipy.stats import spearmanr
import sys, os


def parse_clause(panels_df, idx):
    panel = panels_df['Clause'][idx]
    
    return [x for x in panel.split(',')]


def parse_regional_allele_frequencies(panels_df, idx):
    exclude_cols = ['Generating Anc', 'Clause', 'All Probands']
    regional_freqs = panels_df[panels_df.columns.difference(exclude_cols)].loc[idx]
    
    return regional_freqs.apply(get_allele_frequency)
    
    
def get_allele_frequency(freq_string):
    numerator, denominator = [float(x) for x in freq_string.split('/')]
    
    return numerator / denominator


def get_region_pop_size(panels_df):
    exclude_cols = ['Generating Anc', 'Clause']
    regional_freqs = panels_df[panels_df.columns.difference(exclude_cols)].loc[0]
    pop_sizes = regional_freqs.apply(lambda x: int(x.strip().split('/')[1]))
    
    return pop_sizes


def plot_heatmap(true_freqs, est_freqs, outfiles, title, bins, log_freq_range, max_count):
    ## Calculate Spearman correlation
    spear, p = spearmanr(true_freqs, est_freqs)
    spear = np.round(spear, 3)

    matplotlib.rcParams.update({'font.size': 16})
    epsilon = 0.00001
    true_freqs = np.log10(np.array(true_freqs) + epsilon)
    est_freqs = np.log10(np.array(est_freqs) + epsilon)

    means = []
    bin_centres = []
    for i in range(len(bins)-1):
        print bins[i], bins[i+1]
        ix = np.where((est_freqs >= bins[i]) & (est_freqs < bins[i+1]))[0]
        print ix
        if len(ix) > 0:
            exp_freqs = 10 ** true_freqs[ix]
            print np.log10(np.mean(exp_freqs))
            means.append(np.log10(np.mean(exp_freqs)))
            bin_centres.append((bins[i] + bins[i+1]) / 2.)

    # import IPython; IPython.embed()
    norm=LogNorm(vmin=1, vmax=max_count)
    plt.hist2d(est_freqs, true_freqs, cmap='viridis', bins=bins, norm=norm)
    plt.plot(log_freq_range, log_freq_range, color='lightgrey', linestyle='--')
    plt.plot(bin_centres, means, 'o', color='orange', markersize=5, markerfacecolor="None")
    # plt.hist2d(true_freqs, est_freqs, cmap='viridis', bins=bins, norm=norm)
    # plt.plot(log_freq_range, log_freq_range, color='lightgrey', linestyle='--')
    # plt.plot(means, bin_centres, 'o', color='orange', alpha=0.5)
    plt.xlim(log_freq_range)
    plt.ylim(log_freq_range)
    plt.xticks(np.arange(-5, 0))
    plt.yticks(np.arange(-5, 0))
    plt.xlabel('Estimated allele frequency (log scale)')
    plt.ylabel('True allele frequency (log scale)')
    plt.title(title + '\nSpearman correlation: ' + str(spear))
    cb = plt.colorbar()
    cb.set_label('Number of regions')
    plt.tight_layout()
    plt.axis('equal')

    for outfile in outfiles:
        plt.savefig(outfile)

    plt.clf()



def read_freqs(freq_outfile):
    with tables.open_file(freq_outfile) as f:
        freqs = f.root.regional_expected_freqs[:]

    freq_df = pd.DataFrame(freqs, columns=['Region', 'Estimate', 'N'])
    freq_df['Estimate'] = freq_df['Estimate'].astype(float)
    freq_df['N'] = freq_df['N'].astype(float)

    return freq_df


def load_freq_accuracy(true_freqs_file, est_freqs_file, plot_out_file,
                        est_format):
    true_freqs = pd.read_csv(true_freqs_file)
    if est_format == 'txt':
        est_freqs = pd.read_csv(est_freqs_file)
    else:
        est_freqs = read_freqs(est_freqs_file)
    true_freqs['freq'] = true_freqs['nobs'] / true_freqs['N']
    est_freqs['freq'] = est_freqs['Estimate'] / true_freqs['N']

    total_est = est_freqs[est_freqs['Region'] == 'All Probands']['freq'].values[0]
    total_true = true_freqs[true_freqs['region'] == 'All Probands']['freq'].values[0]
    inflation = total_est / total_true

    corrected_true = true_freqs[true_freqs['region'] != 'All Probands']['freq'].values
    corrected_est = est_freqs[est_freqs['Region'] != 'All Probands']['freq'].values / inflation

    # corrected_true = corrected_true[corrected_true['region'] != 'All Probands'].values
    # corrected_est = corrected_est[corrected_est['Region'] != 'All Probands'].values

    return corrected_true, corrected_est


def plot_freq_accuracy(true_freqs_file, est_freqs_file, plot_out_file,
                        est_format):
    true_freqs = pd.read_csv(true_freqs_file)
    if est_format == 'txt':
        est_freqs = pd.read_csv(est_freqs_file)
    else:
        est_freqs = read_freqs(est_freqs_file)
    true_freqs['freq'] = true_freqs['nobs'] / true_freqs['N']
    est_freqs['freq'] = est_freqs['Estimate'] / true_freqs['N']

    # all_true_counts = true_freqs['nobs'].values
    # all_est_counts = est_freqs['Estimate'].values
    total_est = est_freqs[est_freqs['Region'] == 'All Probands']['freq'].values[0]
    total_true = true_freqs[true_freqs['region'] == 'All Probands']['freq'].values[0]
    inflation = total_est / total_true

    all_true_counts = true_freqs['freq'].values
    all_est_counts = est_freqs['freq'].values / inflation

    plt.plot([0, 1], [0, 1], '--', color='lightgrey')
    plt.scatter(all_true_counts, all_est_counts)
    plt.title('True vs Estimated Regional Allele Frequencies\nInflation Factor: ' +\
              str(np.round(inflation, 2)) + ', Global allele frequency: ' +\
              str(np.round(total_true, 4)))
    plt.xlabel('True Allele Frequency')
    plt.ylabel('Estimated Allele Frequency')
    max_val = np.max([np.max(all_true_counts), np.max(all_est_counts)])
    plt.xlim(0, max_val)
    plt.ylim(0, max_val)
    
    plt.savefig(plot_out_file)
    plt.clf()


def main(args):
    ##--------------------------------------------------
    ## ISGen plots
    ##--------------------------------------------------
    outpath = os.path.expanduser(args.outpath)
    args.true_freqs_files = os.path.expanduser(args.true_freqs_files)
    args.est_freqs_files = os.path.expanduser(args.est_freqs_files)
    true_freqs_files = [line.strip() \
                        for line in open(args.true_freqs_files, 'r')]
    est_freqs_files = [line.strip() \
                        for line in open(args.est_freqs_files, 'r')]
    assert len(true_freqs_files) == len(est_freqs_files)
    plot_out_files = [outpath + '/plot_' + str(i) + '.png' \
                        for i in range(len(est_freqs_files))]

    true_ISGen_vals = []
    est_ISGen_vals = []
    for true, est, out in zip(true_freqs_files, est_freqs_files, plot_out_files):
        print "Loading", out
        new_true_ISGen_vals, new_est_ISGen_vals = load_freq_accuracy(true, est, out, args.est_freqs_format)
        true_ISGen_vals.extend(new_true_ISGen_vals)
        est_ISGen_vals.extend(new_est_ISGen_vals)



    ##--------------------------------------------------
    ## Kinship plots
    ##--------------------------------------------------
    kinshipfile = "/Users/dnelson/project/anc_finder/data/kinship_BALSAC.csv"
    panelsfile = "/Users/dnelson/project/anc_finder/results/validation_panels/BALSAC_hom_panels.txt"
    regionsfile = "/Users/dnelson/project/anc_finder/data/BALSAC_proband_regions.txt"


    kinship = pd.read_csv(kinshipfile, sep=" ", index_col=0)
    panels = pd.read_csv(panelsfile, sep="\t")
    regions = pd.read_csv(regionsfile, sep="\t", header=None, names=["Region"], index_col=0)
    regions['Region'] = pd.core.strings.str_strip(regions['Region'])

    true_kinship_vals = []
    est_kinship_vals = []
    for panel_num in range(100):
        print "Loading panel", panel_num
        panel_kinship = kinship[parse_clause(panels, panel_num)]
        panel_kinship = pd.concat([panel_kinship, regions], axis=1, join_axes=[panel_kinship.index])

        carrier_rates = pd.DataFrame()
        carrier_rates['Estimated'] = panel_kinship.groupby('Region').mean().sum(axis=1) / 2
        carrier_rates['True'] = parse_regional_allele_frequencies(panels, panel_num).to_frame()
        carrier_rates['N'] = get_region_pop_size(panels)
        total_true = (carrier_rates['True'] * carrier_rates['N']).sum()
        total_estimated = (carrier_rates['Estimated'] * carrier_rates['N']).sum()
        inflation_factor = total_estimated / total_true
        carrier_rates['Estimated'] = carrier_rates['Estimated'] / inflation_factor

        true_kinship_vals.extend(carrier_rates['True'].values)
        est_kinship_vals.extend(carrier_rates['Estimated'].values)

    print len(true_kinship_vals), len(est_kinship_vals)


    log_freq_range = [-5, -1]
    bins = np.arange(log_freq_range[0], log_freq_range[1], 0.055)
    max_kinship_count = np.max(np.histogram2d(np.log10(np.array(true_kinship_vals)), np.log10(np.array(est_kinship_vals)), bins=bins)[0])
    max_ISGen_count = np.max(np.histogram2d(np.log10(np.array(true_ISGen_vals)), np.log10(np.array(est_ISGen_vals)), bins=bins)[0])
    max_count = np.max([max_kinship_count, max_ISGen_count])


    outfiles = [
            outpath + 'freq_kinship_hist.png',
            outpath + 'freq_kinship_hist.svg',
            outpath + 'freq_kinship_hist.ps',
            outpath + 'freq_kinship_hist.pdf']
    plot_heatmap(true_kinship_vals, est_kinship_vals, outfiles, 'Kinship', bins, log_freq_range, max_count)

    outfiles = [
            outpath + 'freq_isgen_hist.png',
            outpath + 'freq_isgen_hist.svg',
            outpath + 'freq_isgen_hist.ps',
            outpath + 'freq_isgen_hist.pdf']
    plot_heatmap(true_ISGen_vals, est_ISGen_vals, outfiles, 'ISGen', bins, log_freq_range, max_count)


if __name__ == "__main__":
    args = argparse.Namespace(
            true_freqs_files='~/project/anc_finder/results/validation_briaree/300K/true_freqs_files_tournesol.txt',
            est_freqs_files='~/project/anc_finder/results/validation_briaree/300K/est_freqs_files_tournesol.txt',
            est_freqs_format='hdf5',
            outpath='~/temp/'
            # outpath='~/Documents/papers/allele-sim/resources/'
            )
    main(args)
