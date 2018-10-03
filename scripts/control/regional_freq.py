from __future__ import division
import sys, os
import control_calc as calc
import numpy as np
from functools import reduce
import scipy.stats
import warnings
import argparse
import tables
import scipy.sparse as sparse
from collections import defaultdict
from profilehooks import timecall


def load_regions_h5(contrib_h5file):
    ## Connect to file
    with tables.open_file(contrib_h5file, 'r') as f:
        ## Load all regions to make sure they are the same
        nodes = f.list_nodes(f.root.sparse_hist)
        regions = sorted([node.name.strip() for node in nodes])

    return regions


def load_sparse_contribs(node):
    """ Loads contribs stored in sparse.coo_matrix format in hdf5 node """
    ninds_region = node._v_attrs['ninds']
    r, c, d = np.transpose(np.asarray(node[:])).reshape(3, -1)
    contribs = sparse.csr_matrix((d, (r, c)))

    return contribs, ninds_region


def normalize_sparse_rows(sp_array):
    """
    Normalizes a sparse array along its rows, returning another sparse array.
    """
    assert type(sp_array) is scipy.sparse.coo.coo_matrix

    ## Sum across each row
    totals = np.array(sp_array.sum(axis=1)).ravel()

    ## Divide each entry by the row sum
    new_dat = [1. * d / totals[r] for d, r in zip(sp_array.data, sp_array.row)]

    ## Build new sparse array with the normalized data
    norm = scipy.sparse.coo_matrix((new_dat, (sp_array.row, sp_array.col)))

    return norm


def sparse_expectation(sp_array, weights):
    """
    Calculates the expectation of each row (a distribution) of a sparse array,
    where the support of the distribution is given by the column indices.
    Weights row expectations by values provided in `weights`
    """
    nrows = np.max(sp_array.row) + 1
    row_expectations = np.zeros(nrows)

    for r, c, d in zip(sp_array.row, sp_array.col, sp_array.data):
        row_expectations[r] += c * d * weights[r]

    return row_expectations


def get_genotype_weights(genotypes):
    """
    Returns weights for allele frequency distributions based on the
    probability that each individual would have inherited an allele from their
    parent in the tree.
    """
    ## If tree ind is a het, their offspring have only a 50 percent
    ## chance of receiving an allele. If tree ind is a hom, we don't need
    ## to adjust the distribution, which was simulated with all hets
    ## NOTE: Underestimates for inds with both parents in tree +p3 id:118
    weights = np.ones(len(genotypes))
    for i, genotype in enumerate(genotypes):
        if np.max(genotype) == 1:
            weights[i] = 0.5

    return weights


def get_tree_expectations(T, all_contribs):
    """
    Returns a list of allele frequency distributions for the individuals
    in each tree in T, along with the parent genotypes of each ind in the
    boundary of the trees.
    """
    print "Calculating expectations of trees..."

    ## Store ind contributions in a dict
    ind_contrib_dict = {}
    tree_expectations = np.zeros(len(T.alltrees))

    for i, tree in enumerate(T.alltrees):
        if i % 1000 == 0:
            print i, '/', len(T.alltrees)
        ## Get the boundary inds who have not been calculated already,
        ## and the genotypes of their parents
        boundary_dict = T.getboundary(i)
        boundary_inds, genotypes = zip(*boundary_dict.items())
        new_inds = [ind for ind in boundary_inds if ind not in ind_contrib_dict]

        if len(new_inds) > 0:
            ## Locate boundary inds in indlist
            boundary_indices = np.array([T.ind_dict[ind] for ind in new_inds])

            ## Group ind contribs for the tree and normalize
            raw_contribs = all_contribs[boundary_indices].tocoo()
            contribs = normalize_sparse_rows(raw_contribs)

            ## Calculate expected regional allele frequencies of each tree
            boundary_weights = get_genotype_weights(genotypes)
            expectations = sparse_expectation(contribs, boundary_weights)

            ## Store new expectations in contrib dict
            new_contribs = {ind: e for ind, e in
                                zip(new_inds, expectations)}
            ind_contrib_dict.update(new_contribs)

        ## Sum expectations of boundary inds to get total for the tree
        tree_expectation = [ind_contrib_dict[ind] for ind in boundary_inds]
        tree_expectations[i] = np.sum(tree_expectation)

    return tree_expectations


def get_regional_expectations(T, regions, contrib_h5):
    ## open allele dropping contribution results
    ## TODO: Cleaner if this file isn't open for the whole loop +p4 id:163
    f = tables.open_file(contrib_h5, 'r')
    node_dict = {node.name.strip(): node
                    for node in f.list_nodes(f.root.sparse_hist)}

    ## store freqs and number of probands in each region in a dict
    expected_freqs = {}
    ninds_dict = {}

    for region in regions:
        print "loading contributions to region", region
        ## data stored in sparse format as row, col, dat
        region_node = node_dict[region]
        all_contribs, ninds_region = load_sparse_contribs(region_node)
        ninds_dict[region] = ninds_region

        ## Load regional contributions and parent genotypes of inds in the
        ## boundary of each tree
        expected_freqs[region] = get_tree_expectations(T, all_contribs)

    f.close()

    return expected_freqs, ninds_dict


def get_E_rF(E_F, E_r, freq_params):
    """
    Returns approximation of E[r | \Gamma, F] for tree with global
    expected allele count E_F and regional expected allele count E_r
    """
    if E_F == 0:
        return 0

    assert type(freq_params) is str
    N, n, k = [int(x.strip()) for x in freq_params.split(',')]
    F_obs = 1. * k / n * N

    return E_r * F_obs / E_F


def reg_expectation(T, regions, control_liks, contrib_h5, freq_params):
    """
    Calculates regional frequency expectations based on the pre-calculated
    likelihoods of each tree having produced the observed allele frequency
    """
    regional_expectations = []
    expected_freqs, ninds_dict = get_regional_expectations(T,
                                                           regions,
                                                           contrib_h5)
    global_expectations = expected_freqs['All Probands']

    for region in regions:
        print "Calculating expected allele frequencies for region", region

        regional_freqs = expected_freqs[region]
        numerator = 0.
        denominator = 0.

        for E_r, c, t, E_F in zip(regional_freqs, control_liks, T.climb_liks,
                                global_expectations):
            E_rF = get_E_rF(E_F, E_r, freq_params)
            numerator += E_rF * c * t
            denominator += c * t

        if numerator == 0:
            ## Avoids division by zero if expected frequency is zero
            regional_expectations.append(0)
        else:
            assert denominator != 0
            regional_expectation = numerator / denominator
            regional_expectations.append(regional_expectation)

    return regional_expectations, ninds_dict, expected_freqs


def load_tree_freq_liks(control_freq_likfile):
    with tables.open_file(control_freq_likfile, 'r') as f:
        control_freqs = f.root.control_liks[:]

    return map(lambda x: 2.**x, control_freqs)


def write_regional_freqs(outfile, regions, reg_freqs, tree_freqs,
                                                         ninds_dict):
    """
    Writes expected regional allele frequencies over all simulated trees,
    as well as for each individual tree.
    """
    ## Filters for compression
    filters = tables.Filters(complevel=5, complib='blosc')

    with tables.open_file(outfile, 'w') as f:
        ## Write expected allele frequency for each tree for the region
        for region, freqs in tree_freqs.iteritems():
            ninds = ninds_dict[region]

            ## Suppress warning about naming convention
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ca = f.create_carray(f.root, region, tables.FloatAtom(),
                            shape=tree_freqs[region].shape, filters=filters)
            ca[:] = tree_freqs[region]
            ca.attrs.ninds = ninds

        ## Convert array values to string for writing to character array
        ninds_array = np.asarray([str(ninds_dict[r]) for r in regions])
        reg_freqs = map(str, reg_freqs)
        region_freqs_array = np.asarray(zip(regions, reg_freqs, ninds_array))

        ## Write total expected allele frequency per region
        title = 'region\texpected_alleles\tnum_inds'
        f.create_array(f.root, 'regional_expected_freqs', region_freqs_array,
                                title=title)


def main(args):
    ## Parse CLIs
    pedfile = os.path.expanduser(args.pedfile)
    contrib_h5 = os.path.expanduser(args.contrib_h5)
    climb_likfile = os.path.expanduser(args.climb_likfile)
    outfile = os.path.expanduser(args.outfile)
    control_freq_likfile = os.path.expanduser(args.control_freq_likfile)

    ## Load regions if provided, otherwise calculate expected allele
    ## frequencies for all regions present in the contrib_h5 file
    if args.regionsfile is not None and args.regionsfile != 'None':
        regions = []
        with open(args.regionsfile, 'r') as f:
            for line in f:
                regions.append(line.strip('\r\n'))
        print "Loaded regions:", regions
    else:
        regions = load_regions_h5(contrib_h5)

    T = calc.Tree(pedfile, climb_likfile)

    ## Store the indices of all individuals in the pedigree
    ind_dict = dict(zip(T.inds, np.arange(len(T.inds))))

    print "Loading tree frequency likelihoods..."
    control_liks = load_tree_freq_liks(control_freq_likfile)

    print "Calculating expected regional allele frequencies..."
    reg_freqs, ninds_dict, tree_freqs = reg_expectation(T,
                        regions, control_liks, contrib_h5, args.freq_params)

    ## Write regional frequency estimates to file, both total and per-tree
    print "Writing to file..."
    write_regional_freqs(outfile, regions, reg_freqs, tree_freqs, ninds_dict)

    print "Finished!"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--regionsfile", metavar='|',
                        help="""File listing regions for which to calculate
                                expected allele frequencies""")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-p", "--pedfile", metavar='|',
                        help="Pedigree file",
                        required=True)
    requiredNamed.add_argument("-c", "--contrib-h5", metavar='|',
                        help="Path to hdf5 file containing regional " +\
                            "contribution of inds in pedigree",
                        required=True)
    requiredNamed.add_argument("-l", "--climb_likfile", metavar='|',
                        help="File containing the results of climbing" +\
                        " simulations",
                        required=True)
    requiredNamed.add_argument("-f", "--control-freq-likfile", metavar='|',
                        help="Path to file giving likelihood of each tree" +\
                        " having produced the observed allele frequency " +\
                        "in the entire population",
                        required=True)
    requiredNamed.add_argument("-o", "--outfile", metavar='|',
                        help="File to store simulation output in csv format",
                        required=True)
    requiredNamed.add_argument("-F", "--freq_params", metavar='|',
                        help="Parameters describing observations of the " +\
                            "population allele frequency, in the format " +\
                            "'pop_size,num_sample,num_obs_carriers'",
                        required=True)

    args = parser.parse_args()

    main(args)
