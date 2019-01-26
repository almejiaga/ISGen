import sys, os
sys.path.append(os.path.abspath('../climb/'))
sys.path.append(os.path.abspath('../'))
import numpy as np
from collections import Counter, defaultdict
import argparse
import copy
import pytest
import tempfile
import ped, climb, run


## Default args to be used in test simulations
args = argparse.Namespace(label="test",
                iterations=1000,
                sim_homs=0.5,
                hetfile='test_data/hets2.txt',
                homfile='test_data/homs2.txt',
                proband_regions='test_data/regions2.txt',
                contrib_outfile=tempfile.NamedTemporaryFile().name,
                climb_outfile=tempfile.NamedTemporaryFile().name,
                control_outfile=tempfile.NamedTemporaryFile().name,
                regional_freq_outfile=tempfile.NamedTemporaryFile().name,
                contrib_iterations=1000,
                delimiter=None,
                compressed=False,
                verbose=False,
                climb_iterations=1000,
                freq_params='100,10,1',
                regions_to_calculate=None)

## Specific cases to test
params = []
kinship_vals = [True, False]

pedfile2 = 'test_data/pedEx2.txt'
genotypes2 = defaultdict(int, {1: 1, 2: 1, 3: 1})
sim_homs_vals = [0, 0.5, 1]
for sim_homs in sim_homs_vals:
    for kinship in kinship_vals:
        params += [(pedfile2, genotypes2, sim_homs, kinship)]

pedfile4 = 'test_data/pedEx4.txt'
genotypes4 = defaultdict(int, {1: 1, 2: 1, 3: 1, 4: 1})
sim_homs_vals = [0, 1]
for sim_homs in sim_homs_vals:
    for kinship in kinship_vals:
        params += [(pedfile4, genotypes4, sim_homs, kinship)]



## Define scope as entire module, so the results can be reused by all tests
@pytest.fixture(scope='module', params=params)
def sim_data(request):
    """ Generate simulation data to be used in other tests """
    ## Set parameters which vary between runs
    args.pedfile, known_genotypes, args.sim_homs, args.kinship = request.param

    ## Set up simulations
    P = ped.Pedigree(args.pedfile)
    ind_cones = P.allowedinds(known_genotypes.keys())[2]
    parent_dict = P.parent_dict.copy()

    signals = dict()
    signals[0] = 0
    for ind in parent_dict:
        signals[ind] = 0
                    
    ancs_dict = P.ancestors_dict(known_genotypes.keys())

    A = climb.AncFinder(ind_cones, parent_dict, ancs_dict, args)

    ## Run Simulations
    ancs, liks, trees = [], [], []
    for i in range(args.iterations):
        genotypes = known_genotypes.copy()
        sig = signals.copy()

        ## Simulate coalescence of all affected alleles
        anc, lik, tree = A.coalesce(genotypes, sig)
        ancs.append(anc)
        liks.append(lik)
        trees.append(tree)

    yield {'ancs': ancs, 'liks': liks, 'args': args, 'trees': trees}


def test_anc_symmety(sim_data):
    """
    Two members of a couple, who have offspring only with each other, should
    be chosen equally often. Tests for equality using 95 percent confidence
    intervals for binomial distribution.
    """
    ## Don't need to run over all parameter values
    if sim_data['args'].sim_homs == 0.5:
        counts = Counter(sim_data['ancs'])

        binomial_sd = np.sqrt(sim_data['args'].iterations * 0.25)
        mean = np.mean(counts.values())
        conf95 = (mean - 2 * binomial_sd, mean + 2 * binomial_sd)

        for count in counts.values():
            assert conf95[0] < count < conf95[1]


def test_homs(sim_data):
    """ Tests whether homozygotes can be avoided or forced when possible """
    ## Make sure no homozygotes are created if we specify probability 0
    if sim_data['args'].sim_homs == 0:
        for tree in sim_data['trees']:
            assert 2 not in set(tree.values())

    ## Make sure all possible homozygotes are created if we specify
    ## probability 1
    elif sim_data['args'].sim_homs == 1:
        for tree in sim_data['trees']:
            if sim_data['args'].pedfile == 'test_data/pedEx2.txt':
                assert tree[201] == 2 or tree[202] == 2
            elif sim_data['args'].pedfile == 'test_data/pedEx4.txt':
                assert tree[6] == 2 and tree[7] == 2
                assert tree[9] == 1 and tree[10] == 1
                assert tree[11] == 2 or tree[12] == 2


def test_liks(sim_data):
    """ Checks likelihoods which depend on the validity of homozygotes """
    for lik, tree in zip(sim_data['liks'], sim_data['trees']):
        if sim_data['args'].pedfile == 'test_data/pedEx2.txt':
            if sim_data['args'].sim_homs == 0:
                assert set(sim_data['liks']) == set([-8])
            elif sim_data['args'].sim_homs == 0.5:
                assert set(sim_data['liks']) == set([-7])
            elif sim_data['args'].sim_homs == 1:
                assert set(sim_data['liks']) == set([-8])
        elif sim_data['args'].pedfile == 'test_data/pedEx4.txt':
            if sim_data['args'].kinship is True:
                ## TODO: NotImplemented
                pass
            elif sim_data['args'].sim_homs == 0:
                if tree[9] == 0 or tree[10] == 0:
                    # assert lik = -7 # if no parent weights used
                    assert np.isclose(lik, -7 + np.log2(0.5 / 0.6))
                else:
                    # assert lik = -8 # if no parent weights used
                    assert np.isclose(lik, -8 + np.log2(0.5 / 0.4))
            elif sim_data['args'].sim_homs == 1:
                assert set(sim_data['liks']) == set([-10])


def test_CLI(sim_data):
    """ Validates how args are passed without actually using the CLI """
    run.main(sim_data['args'])
