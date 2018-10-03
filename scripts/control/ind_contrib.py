# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 21:14:44 2015

@author: dom
"""
from __future__ import division
import sys,os
import numpy as np
from collections import defaultdict
import argparse
import csv
import scipy.sparse
import itertools
import tables
import warnings
import control_calc as calc
import hdf5_to_sparse


def initialize_hf_arrays(tot_ninds, region_dict, niter, h5file):
    ##TODO Might want to check if file exists already... will overwrite
    f = tables.open_file(h5file, 'w')
    g = f.create_group(f.root, 'raw')
    fill = np.zeros((tot_ninds, niter))

    ## Create the blank contrib file
    for region, inds in region_dict.iteritems():
        filters = tables.Filters(complevel=5, complib='blosc')
        ## Suppress warning about naming convention
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ca = f.create_carray(g, region,
                                tables.Atom.from_dtype(fill.dtype),
                                shape=(fill.shape), filters=filters)
        ca[:] = fill
        ca.attrs.ninds = len(inds)

    ## Close the file
    f.close()

    return None


def increment_h5_counts(regions, iteration, sparse_contribs, h5file):
    f = tables.open_file(h5file, 'a')
    sparse_contribs = sparse_contribs.tocoo()

    try:
        ## Append the new contributions
        for i, j, k, in itertools.izip(sparse_contribs.row, sparse_contribs.col,
                                                        sparse_contribs.data):
            rg = f.get_node(f.root.raw, regions[j])
            rg[i, iteration] = k
    except:
        ## Make sure the file is closed if anything goes wrong
        f.close()
        raise

    ## Close the file
    f.close()

    return None


def load_regions(regionfile, delimiter=None):
    import csv

    if delimiter is None:
        delimiter = '\t'

    region_dict = defaultdict(list)

    with open(regionfile, 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for line in reader:
            try:
                region_dict[line[1]].append(int(float(line[0])))
            except IndexError:
                print("Error parsing file: try using --delimiter option")
                sys.exit()

    return region_dict


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


def main(args):
    pedfile = os.path.expanduser(args.pedfile)
    outfile = os.path.expanduser(args.outfile)
    iterations = args.iterations

    ## Set verbose flag
    vprint.verbose = args.verbose

    vprint(("Loading pedigree..."))
    P = calc.Pedigree(pedfile)
    path = P.buildpath()
    vprint(("Dropping alleles from all individuals..."))

    if args.regionfile is not None:
        regionfile = os.path.expanduser(args.regionfile)
        vprint(("Loading proband regions..."))
        region_dict = load_regions(regionfile)
        vprint(("Done"))
        region_dict['All Probands'] = P.probands
    else:
        regionfile = None
        region_dict = {'All Probands': P.probands}

    regions = sorted(region_dict.keys())

    ## Create the database which will be written to every iteration.
    initialize_hf_arrays(len(P.inds), region_dict, iterations, outfile)

    for k in range(iterations):
        vprint(k + 1, " / ", iterations)
        vprint(("Passing alleles through pedigree..."))
        alleles = P.passalleles(path)
        vprint(("Done."))

        vprint(("Collecting allele counts per region..."))
        contribs = P.get_contrib(alleles, region_dict)
        ## Regions are sorted so that indices match with values returned from
        ## simulations in control_calc.py
        vprint(("Done."))

        vprint(("Adding allele counts to total"))
        ## Increase each (ind, region, value) entry by 1.
        increment_h5_counts(regions, k, contribs, outfile)
        vprint(("Done."))

    vprint(("Allele dropping complete."))

    ## Write sparse allele frequency histograms to same output file
    sparse_args = argparse.Namespace(dropfile=outfile, verbose=args.verbose)
    hdf5_to_sparse.main(sparse_args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--regionfile", metavar='',
                        help="File listing the region of each proband, " +\
                            "with no header")
    parser.add_argument("-m", "--delimiter", metavar='',
                        help="Delimiter used in --regionfile. If space, " +\
                        "surround by double quotes. Default is tab.")
    parser.add_argument("-v", "--verbose", action="store_true")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-p", "--pedfile", metavar='',
                        help="Pedigree file",
                        required=True)
    requiredNamed.add_argument("-o", "--outfile", metavar='',
                        help="File to store simulation output in hdf5 format",
                        required=True)
    requiredNamed.add_argument("-i", "--iterations", metavar='',
                        help="Number of allele drops to perform from each " +\
                            "individual in the pedigree",
                        type=int, required=True)

    args = parser.parse_args()


    main(args)
