import sys, os
sys.path.append(os.path.abspath('../control/'))
sys.path.append(os.path.abspath('../climb/'))
import numpy as np
from collections import Counter, defaultdict
import argparse
import copy
import pytest
import tempfile
import tables
import scipy.sparse
import ind_contrib, ped


## Default args to be used in test simulations
args = argparse.Namespace(regionfile="test_data/regions2.txt",
            pedfile='test_data/pedEx2.txt',
            delimiter=None,
            compressed=None,
            verbose=False,
            iterations=1000)

## Specific cases to test - with and without regions defined
regionfiles = ["test_data/regions2.txt", None]
outfiles = [tempfile.NamedTemporaryFile().name for r in regionfiles]

params = zip(regionfiles, outfiles)


## Define scope as entire module, so the results can be reused by all tests
@pytest.fixture(scope="module", params=params)
def contrib_data(request):
    """ Generate simulation data to be used in other tests """
    ## Set parameters which vary between runs
    args.regionfile, args.outfile = request.param

    ## Perform allele dropping simulations
    ind_contrib.main(args)

    ## Load results back for testing
    with tables.open_file(args.outfile, 'r') as f:
        ## Get nodes and their names, which correspond to regions in the
        ## probands
        node_dict = {n.name: n for n in f.list_nodes(f.root.sparse_hist)}
        contrib_dict = {}

        ## Load data from nodes as a sparse array
        for region, node in node_dict.iteritems():
            ## Data is stored in three-row sparse format, as row, col, dat
            row, col, dat = np.transpose(np.asarray(node[:])).reshape(3, -1)

            ## Load indices of each individual in the pedigree, so accessing
            ## their contributions later can be done by their label
            ped_indices = np.genfromtxt(args.pedfile, usecols=(0), dtype=int,
                                        skip_header=1)
            new_rows = [ped_indices[i] for i in row]
            sparse_contribs = scipy.sparse.coo_matrix((dat, (new_rows, col)))

            ## Store contributions by region in a dict
            contrib_dict[region] = sparse_contribs.tocsr()

    yield contrib_dict


def test_counts(contrib_data):
    """
    Checks that allele counts are within 95 percent confidence intervals.
    Data are stored in arrays with format [ind, num_alleles] = count
    """
    ## Find the probands - they should always contribute their own alleles
    P = ped.Pedigree(args.pedfile)

    ## TODO: Test individuals with multiple offspring +p2 id:112
    ## Case when no regions are specified - only all probands
    regions = contrib_data.keys()
    if len(regions) == 1:
        all_contribs = contrib_data[regions[0]]
        for prob in P.probands:
            ## Probands should always contribute 1 allele
            assert all_contribs[prob, 1] == args.iterations

        ## Test specific individuals
        counts = Counter(all_contribs[11])
        binomial_sd = np.sqrt(args.iterations * 0.25)
        mean = np.mean(counts.values())
        conf95 = (mean - 2 * binomial_sd, mean + 2 * binomial_sd)

        for count in counts.values():
            assert conf95[0] < count < conf95[1]


def test_output_exists(contrib_data):
    assert os.path.exists(args.outfile)
