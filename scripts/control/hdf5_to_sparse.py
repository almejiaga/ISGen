from __future__ import division
import scipy.sparse as sparse
import numpy as np
import multiprocessing
import tables
import sys, os
import warnings
import argparse


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
    ## Parse CLI args
    h5 = args.dropfile

    ## Set verbose flag
    vprint.verbose = args.verbose

    h5f = tables.open_file(h5, 'a')
    nodes = h5f.list_nodes(h5f.root.raw)
    ninds = nodes[0].shape[0]

    ## Group for storing sparse output
    sparse_name = 'Histogram of allele contributions per ind by region'
    g = h5f.create_group(h5f.root, 'sparse_hist', sparse_name)

    ## Rows and data for the histogram arrays are constant
    r = range(ninds)
    d = list(np.ones(ninds, dtype=np.int16))

    ## Add data from hdf5 file
    for node in nodes:
        region = node.name
        vprint("Loading data from region", region)
        sp_array = sparse.coo_matrix((1,1))
        nsims = node.shape[1]

        ## Load data from file in chunks to limit memory usage
        blocksize = 5
        sim_slices = range(0, nsims, blocksize)

        r_slice = sorted(r * blocksize)
        d_slice = d * blocksize

        for i in range(1, len(sim_slices)):
            vprint("Loading up to simulation", sim_slices[i])
            loaded_dat = node[:, sim_slices[i-1] : sim_slices[i]]
            vprint("Block shape:", loaded_dat.shape)

            ## Construct new rows, cols, and data for coo_matrix. Each
            ## simulation gets added as a block to the end.
            rows = np.append(sp_array.row, r_slice)
            data = np.append(sp_array.data, d_slice)
            cols = np.append(sp_array.col, loaded_dat.ravel())

            ## Construct sparse coo_matrix from new data.
            sp_array = sparse.coo_matrix((data, (rows, cols)))
            sp_array = sp_array.tocsc().tocoo()
            vprint("Condensed shape:", sp_array.shape)
            vprint(len(sp_array.data), "entries")

        ## Load leftover columns
        vprint("Loading up to simulation", node.shape[1] - 1)
        loaded_dat = node[:, sim_slices[-1] :]
        vprint("Block shape:", loaded_dat.shape)
        num_loaded_sims = loaded_dat.shape[1]
        r_slice = sorted(r * num_loaded_sims)
        d_slice = d * num_loaded_sims

        ## Construct new rows, cols, and data for coo_matrix. Each simulation
        ## gets added as a block to the end.
        rows = np.append(sp_array.row, r_slice)
        data = np.append(sp_array.data, d_slice)
        cols = np.append(sp_array.col, loaded_dat.ravel())

        ## Construct sparse coo_matrix from new data.
        sp_array = sparse.coo_matrix((data, (rows, cols)))
        sp_array = sp_array.tocsc().tocoo()
        vprint("Condensed shape:", sp_array.shape)
        vprint(len(sp_array.data), "entries")

        ## Store data in a variable-length array, mimicking coo_matrix format
        ## as an array of shape (3, :) with rows (row, col, data)
        vprint("Writing output to", h5)
        atom = tables.Int32Atom(shape=(3, len(sp_array.data)))

        ## Suppress warning about naming convention
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            ## Transpose data so max row size is not exceeded
            data = np.vstack([sp_array.row, sp_array.col, sp_array.data])
            data = np.transpose(data)

            ## Store in compressed array
            filters = tables.Filters(complevel=5, complib='blosc')
            ca = h5f.create_carray(g, str(region), tables.IntAtom(),
                                    shape=data.shape, filters=filters)
            ca[:] = data

            ## Store number of inds in region
            raw_node = h5f.get_node(h5f.root.raw, str(region))
            ninds_region = raw_node._v_attrs['ninds']
            h5f.set_node_attr(ca, 'ninds', ninds_region)

        h5f.flush()

    ## Close remaining open files
    h5f.close()
    vprint("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", metavar='', action="store_true")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-d", "--dropfile", metavar='|',
                        help="File containing raw allele dropping results" +\
                        " in hdf5 format",
                        required=True)
    requiredNamed.add_argument("-o", "--outfile", metavar='',
                        help="File to store sparse histograms of allele " +\
                        "counts in hdf5 format",
                        required=True)

    main(args)
