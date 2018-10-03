import sys, os
import numpy as np
import tables
import scipy.sparse
import argparse
from collections import defaultdict


def append_sparse(sparse_array, new_row, new_col, new_dat):
    """ Adds entries to existing sparse array, summing matching indices """
    r = np.hstack([sparse_array.row, new_row])
    c = np.hstack([sparse_array.col, new_col])
    d = np.hstack([sparse_array.data, new_dat])

    return scipy.sparse.csr_matrix((d, (r,c))).tocoo()


def initialize_sparse_array():
    """ Initializes a sparse array with a single 0 entry """
    row = np.zeros(1)
    col = np.zeros(1)
    dat = np.zeros(1)

    return scipy.sparse.coo_matrix((dat, (row, col)))


def combine_partial_contribs(partial_files, outfile):
    """
    Combines the sparse results of several allele dropping simulations
    into a single sparse array.
    """
    ## Dict of sparse contribs per region
    region_contribs = defaultdict(initialize_sparse_array)

    for pf in partial_files:
        print "Loading contribs from", pf
        with tables.open_file(pf, 'r') as f:
            ## Iterate through each region node in the sparse file
            for partial_node in f.list_nodes(f.root.sparse_hist):
                ## Load sparse data
                region = partial_node.name
                r, c, d = np.transpose(partial_node[:])

                ## Add results to total
                new_tot = append_sparse(region_contribs[region], r, c, d)
                region_contribs[region] = new_tot

    ## Write results to file
    with tables.open_file(outfile, 'w') as f:
        filters = tables.Filters(complevel=5, complib='blosc')
        g = f.create_group(f.root, 'sparse_hist')

        for region, contribs in region_contribs.iteritems():
            print "Writing contribs for region:", region
            ## Transpose data so max row size is not exceeded
            r, c, d = contribs.row, contribs.col, contribs.data
            contribs_array = np.transpose(np.vstack([r, c, d]))

            ## Store in compressed array
            ca = f.create_carray(g, region.strip(), tables.IntAtom(),
                                shape=(contribs_array.shape), filters=filters)
            ca[:] = contribs_array

            ## Store number of inds in region
            with tables.open_file(partial_files[0], 'r') as pf:
                raw_node = pf.get_node(pf.root.raw, str(region))
                ninds_region = raw_node._v_attrs['ninds']
            f.set_node_attr(ca, 'ninds', ninds_region)


def main(args):
    ## Load paths of partial output files
    with open(args.outpaths_file, 'r') as f:
        partial_files = [os.path.expanduser(l.strip('\n')) for l in f]

    combine_partial_contribs(partial_files, args.outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--outpaths_file', metavar='|',
                    help='File containing paths of all partial output files',
                    required=True)
    parser.add_argument('-o', '--outfile', metavar='|',
                    help='File to write combined allele dropping output',
                    required=True)

    args = parser.parse_args()

    main(args)
