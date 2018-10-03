# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 00:34:19 2015

@author: dom
"""
from __future__ import division
import sys,os
import numpy as np
import control_calc as calc
import conv_utils
import csv
import scipy.stats
import tables
import argparse
import pandas as pd
from profilehooks import profile, timecall

class VPrint(object):
    def __init__(self):
        """ Easily turn print output on/off """
        self.verbose = False

    def __call__(self, *args):
        if self.verbose:
            # Print each argument separately so caller doesn't need to
            # stuff everything to be printed into a single string
            for arg in args:
                print arg,

            print

## Instantiate printing class globally
vprint = VPrint()


def init_array(outfile, array_name):
    """ Initialize extendable arrays to store control likelihoods """
    with tables.open_file(outfile, 'w') as f:
        title = 'Likelihood of tree having produced the observed allele freq'
        f.create_earray(f.root, array_name, atom=tables.FloatAtom(),
                                title=title, shape=(0,))


def write_posteriors(outfile, posteriors, tot_liks):
    """ Writes posterior probabilities and total liks to hdf5 file """
    with tables.open_file(outfile, 'a') as f:
        title = 'Posterior probabilities of each ancestor assuming a ' +\
                'flat prior\nAnc\tProb\tLik'
        f.create_array(f.root, 'posteriors', posteriors, title=title)

        title = 'Total tree likelihoods'
        f.create_array(f.root, 'tot_liks', tot_liks, title=title)


def load_contribs(contribfile):
    ## Load allele frequency distributions for all individuals, over
    ## all probands
    with tables.open_file(contribfile, 'r') as f:
        global_contrib_node = f.get_node(f.root.sparse_hist, 'All Probands')
        x = np.transpose(global_contrib_node[:]).reshape(3, -1)
        all_hists = scipy.sparse.csr_matrix((x[2], (x[0], x[1])))

    return all_hists


def get_boundary(T, treenum, all_contribs):
    """
    Returns the allele frequency distributions of individuals in the boundary
    of the tree, along with the genotypes of the parents (who are in the
    tree) of these individuals.
    """
    boundary_dict = T.getboundary(treenum)
    inds, genotypes = zip(*boundary_dict.items())

    return inds, genotypes


def get_tot_posteriors(climb_outfile, control_outfile):
    """
    Combines climbing and control likelihoods, and outputs total likelihoods
    along with posterior probabilities of each ancestor.
    """
    with tables.open_file(climb_outfile, 'r') as f:
        ancs = f.root.ancs[:].reshape(-1, 1)
        climb_liks = f.root.liks[:].reshape(-1, 1)

    with tables.open_file(control_outfile, 'r') as f:
        control_liks = f.root.control_liks[:].reshape(-1, 1)
        try:
            hap_log2lik = f.root.haplotype_liks[:].reshape(-1, 1)
            print "Found haplotype likelihoods"
        except tables.NoSuchNodeError:
            print "No haplotype likelihoods found"
            hap_log2lik = np.zeros(control_liks.shape)

    ## Load likelihood data into a pandas dataframe
    lik_data = pd.DataFrame(
            np.hstack([ancs, climb_liks, control_liks, hap_log2lik]),
            columns=['anc', 'climb_log2lik', 'control_log2lik', 'hap_log2lik'])

    ## Combine climbing and control likelihoods, and convert from log2 format
    lik_data['tot_log2lik'] = lik_data['climb_log2lik'] + \
                              lik_data['control_log2lik'] +\
                              lik_data['hap_log2lik']
    lik_data['tot_lik'] = lik_data['tot_log2lik'].apply(np.exp2)

    ## Calculate posteriors and convert to np array for writing to h5 file
    lik_data['tot_lik'] = lik_data['tot_lik'] / len(lik_data)
    posteriors = lik_data[['anc', 'tot_lik']].groupby('anc').sum()
    posteriors['anc'] = posteriors.index
    posteriors['prob'] = posteriors['tot_lik'] / posteriors['tot_lik'].sum()

    ## Change back to log2 likelihood to output alongside the posteriors
    posteriors['tot_log2lik'] = posteriors['tot_lik'].apply(np.log2)
    del posteriors['tot_lik']

    return posteriors.values, lik_data['tot_log2lik'].values


# @profile
def main(args):
    ## Parse CLI args
    pedfile = os.path.expanduser(args.pedfile)
    climb_outfile = os.path.expanduser(args.climb_likfile)
    contribfile = os.path.expanduser(args.contrib_file)
    outfile = os.path.expanduser(args.outfile)
    N, n, k = map(int, args.freq_params.split(','))

    ## Set verbose flag
    vprint.verbose = args.verbose

    vprint("Loading pedigree...")
    T = calc.Tree(pedfile, climb_outfile)
    vprint("Done.")

    ## Load contributions of all inds in pedigree
    ## IDEA: Could calculate likelihood of observed regional freqs +p5 id:115
    all_hists = load_contribs(contribfile)

    ## Initialize extendable arrays to store control likelihoods
    init_array(outfile, array_name='control_liks')

    ## Create a dictionary of indices for all individuals, so we can locate
    ## them in the contrib file
    ind_dict = dict(zip(T.inds, np.arange(len(T.inds))))

    ## Instantiate class for handling storage and retrieval of individual
    ## allele contribution probabilities
    C = conv_utils.IndContribs(all_hists, ind_dict, N, n, k)

    vprint("Starting control likelihood calculations...")
    total_prob = 0
    for tree_num in range(len(T.alltrees)):
        if tree_num % 10 == 0:
            vprint("Calculating control likelihood for tree number", tree_num)

        ## We load the allele frequency distributions of the
        ## boundary individuals only, as they are needed.
        inds, parent_genotypes = get_boundary(T, tree_num, all_hists)

        ## Get allele contribution probabilities for each ind in the boundary
        boundary_control_liks = C.get_allele_contribs(inds, parent_genotypes)

        ## Calculate likelihood of this tree having produced the observed
        ## allele frequency
        tree_control_lik = conv_utils.sum_control_liks(boundary_control_liks, k)

        ## Get the likelihood of the tree having produced k observations out of
        ## a sample of size n
        conv_prob = tree_control_lik[k]

        ## Write control likelihoods to file in log base 2
        ## TODO: Write in blocks to reduce open/close cycles of file +p2  id:177
        with tables.open_file(outfile, 'a') as f:
            f.root.control_liks.append([np.log2(conv_prob)])

    posteriors, tot_liks = get_tot_posteriors(climb_outfile, outfile)
    write_posteriors(outfile, posteriors, tot_liks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", action="store_true")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-p", "--pedfile", metavar='',
                        help="Pedigree file",
                        required=True)
    requiredNamed.add_argument("-l", "--climb-likfile", metavar='',
                        help="File containing results of climbing simulations",
                        required=True)
    requiredNamed.add_argument("-o", "--outfile", metavar='',
                        help="File to store control likelihoods for each " +\
                        "tree, in hdf5 format", required=True)
    requiredNamed.add_argument("-c", "--contrib-file", metavar='',
                        help="File containing results of allele dropping " +\
                            "simulations, for all individuals in --pedfile",
                            required=True)
    requiredNamed.add_argument("-f", "--freq-params", metavar='',
                        help="Parameters describing observations of the " +\
                            "population allele frequency, in the format " +\
                            "'pop_size,num_sample,num_obs_carriers'",
                            required=True)

    args = parser.parse_args()

    main(args)
