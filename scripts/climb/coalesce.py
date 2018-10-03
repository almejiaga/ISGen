import sys, os
import ped, climb
import datetime
import subprocess
from collections import defaultdict
import argparse
import numpy as np
import pandas as pd
import tables
import matplotlib.pyplot as plt
import time
import copy
import random

from collections import OrderedDict


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

            print('')

## Instantiate printing class globally
vprint = VPrint()


def make_if_not_exists(path):
    """
    Creates 'path' unless it exists, in which case do nothing.
    """
    try:
        os.makedirs(path)
    except OSError as e:
        ## If error raised because path exists, do nothing
        if e.errno != 17:
            raise


def incremental_write(ancs, liks, trees, blocksize, outfile):
    """ Append to existing csv written from Pandas dataframe """
    simdat = {'Anc': ancs[-blocksize:], 'Log2Lik': liks[-blocksize:],
                    'Tree': trees[-blocksize:]}

    ## Write all output to hdf5 file, in separate nodes
    with tables.open_file(outfile, 'a') as f:
        f.root.ancs.append(simdat['Anc'])
        f.root.liks.append(simdat['Log2Lik'])

        ## Trees, which are variable-length, must be added individually
        for tree in simdat['Tree']:
            tree_inds = tree.keys()
            genotypes = tree.values()
            f.root.trees.append(tree_inds)
            f.root.genotypes.append(genotypes)


def init_output(outfile):
    """
    Creates a timestamped output directory and writes header to output file.
    """
    with tables.open_file(outfile, 'w') as f:
        ## Create extendable arrays so we can incrementally write output
        f.create_earray(f.root, 'ancs', atom=tables.IntAtom(), shape=(0,))
        f.create_earray(f.root, 'liks', atom=tables.FloatAtom(), shape=(0,))

        ## Trees, which are variable-length, must be added individually
        f.create_vlarray(f.root, 'trees', atom=tables.IntAtom())
        f.create_vlarray(f.root, 'genotypes', atom=tables.IntAtom())


def load_genotypes(args):
    """"
    Construct dict of known alleles to be climbed, from files containing
    homs and hets
    """
    known_genotypes = defaultdict(int)
    if args.homfile is not None and args.homfile != "None":
        homfile = os.path.expanduser(args.homfile)
        homs = map(int,[line.strip() for line in open(homfile, 'rU')])
        for hom in homs:
            known_genotypes[hom] = 2

    if args.hetfile is not None and args.hetfile != "None":
        hetfile = os.path.expanduser(args.hetfile)
        hets = map(int,[line.strip() for line in open(hetfile, 'rU')])
        for het in hets:
            known_genotypes[het] = 1

    return known_genotypes


def get_posteriors(simdat):
    """ Calculate posterior probability based on a flat prior """
    simdat['Lik'] = simdat['Log2Lik'].apply(np.exp2)
    ##TODO: Should this be used? +t1
    normalization_factor = simdat['Lik'].mean()
    simdat['Lik'] = simdat['Lik'] / len(simdat)
    posteriors = simdat[['Anc', 'Lik']].groupby('Anc').sum()
    posteriors['Anc'] = posteriors.index
    posteriors['Prob'] = posteriors['Lik'] / posteriors['Lik'].sum()

    ## Change back to log2 likelihood to output alongside the posteriors
    posteriors['Log2Lik'] = posteriors['Lik'].apply(np.log2)
    del posteriors['Lik']

    return posteriors.values


def simulate(args, initial_genotypes, outfile):
    """
    Performs the specified number of climbing simulations, returning
    ancestors, importance sampling likelihoods, and inheritance histories.
    """
    ## Load the pedigree and parent-offspring relationships
    vprint("Loading pedigree...")
    P = ped.Pedigree(os.path.expanduser(args.pedfile))
    
    parent_dict = P.parent_dict
    ##TODO Pass the whole dictionary to improve importance sampling +t1
    ind_cones = P.allowedinds(initial_genotypes.keys())[2]
    
    signals = dict()
    signals[0] = 0
    
    print('create signal dict')
    for ind in parent_dict:
        signals[ind] = 0
                    
    print('create ancestors dict')
    ancs_dict = P.ancestors_dict(initial_genotypes.keys())
    
    print('initialisation')
    A = climb.AncFinder(ind_cones, parent_dict, ancs_dict, args)
    
    ## Perform climbing simulations
    vprint("Performing simulations...")
    blocksize = 100
    simdat_dict = {'Anc':[],'Log2Lik':[], 'Tree':[]}
    for i in range(args.iterations):
        if i % 100 == 0:
            vprint(i, "/", args.iterations)

        ## Make a copy of the initial genotypes so they don't get modified
        genotypes = initial_genotypes.copy()
        sig = signals.copy()

        ## Simulate coalescence of all affected alleles
        anc, lik, tree = A.coalesce(genotypes, sig)
        simdat_dict['Anc'].append(anc)
        simdat_dict['Log2Lik'].append(lik)
        simdat_dict['Tree'].append(tree)
        
        ## Write results to output in blocks
        if i % blocksize == 0:
            incremental_write(
                    simdat_dict['Anc'],
                    simdat_dict['Log2Lik'],
                    simdat_dict['Tree'],
                    blocksize,
                    outfile)
        
    ## Write leftover results
    num_remaining = i % blocksize
    incremental_write(simdat_dict['Anc'],
            simdat_dict['Log2Lik'],
            simdat_dict['Tree'],
            num_remaining,
            outfile)

    simdat_i = pd.DataFrame(simdat_dict)
    simdat_i.index.name = 'Sim'
    
    return simdat_i
    

def main(args):
    ## Make sure some probands are specified as either homs or hets
    homhet_err = 'When climbing alleles, at least one of --homfile or ' +\
                                                        '--hetfile is required'
    assert (args.homfile or args.hetfile), homhet_err

    ## Set verbose flag
    vprint.verbose = args.verbose

    ## Initialize output file and store git details if requested
    outfile = os.path.expanduser(args.outfile)
    init_output(outfile)

    ## Load alleles to climb from files
    known_genotypes = load_genotypes(args)

    ## Perform climbing simulations
    simdat = simulate(args, known_genotypes, outfile)

    ## Write raw results and ancestor likelihoods to file, cleaning up partial
    ## output files
    vprint("Writing output...")
    posteriors = get_posteriors(simdat)
    with tables.open_file(outfile, 'a') as f:
        f.create_array(f.root, 'posteriors', posteriors)
    vprint("Done!")
    

if __name__ == "__main__":
    ## TODO: Allow stopping at first coalescence point +p3 id:6
    ## TODO: Allow climbing of unaffected alleles (Simon's idea?) +p3 id:4 gh:2

    ## Parse CLI arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sim-homs", metavar='|', default=0.5, type=float,
                        help="Frequency with which hets become homs when" +\
                        "climbed to. (default: 0.5)")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-k", "--kinship",
                        action='store_true',
                        help="Flag for use of kinship-based importance sampling,\n" +\
                             " which trades per-iteration computation time for\n" +\
                             " faster per-iteration convergence. Performance \n" +\
                             "increase only seen downstream, for control \n" +\
                             "likelihoods or regional allele frequency estimates")

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-p", "--pedfile", metavar='|',
                        help="Pedigree file",
                        required=True)
    
    requiredNamed.add_argument("-o", "--outfile", metavar='|',
                        help="Directory in which to store output",
                        required=True)
    requiredNamed.add_argument("-i", "--iterations", metavar='|',
                        help="Number of allele drops to perform from each " +\
                            "individual in the pedigree",
                        type=int, required=True)
    

    require_one = parser.add_argument_group('require at least one of')
    require_one.add_argument("-H", "--hetfile", metavar='',
                        help="File containing heterozygous probands")
    require_one.add_argument("-O", "--homfile", metavar='',
                        help="File containing homozygous probands")

    args = parser.parse_args()

    main(args)
