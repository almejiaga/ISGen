from __future__ import division
import sys, os
import numpy as np
import scipy.stats
import scipy.sparse
from itertools import izip
import attr
from profilehooks import profile, timecall


@attr.s
class IndContribs(object):
    """
    Explain how dists will be used later
    """
    all_contrib_hists = attr.ib()
    ind_dict = attr.ib()
    N = attr.ib()
    n = attr.ib()
    k = attr.ib()

    known_allele_contribs = attr.ib(default=attr.Factory(dict))

    def get_allele_contribs(self, inds, parent_genotypes):
        """
        Calculates the probability of each ind in `inds` having produced the
        observed `x` alleles in a sample of size `n` in a population of size
        `N`, where x = range(0, k+1). Genotypes in `parent_genotypes` affect
        how likely these inds are to carry an affected allele
        """
        ## Keys are stored as (ind, parent_genotype) tuples
        ind_genotype_pairs = zip(inds, parent_genotypes)

        ## Pull out inds who we have already calculated
        known_pairs = set(self.known_allele_contribs.keys())
        unknown_pairs = list(set(ind_genotype_pairs).difference(known_pairs))

        ## If we have calculated all inds already, we're done
        if len(unknown_pairs) == 0:
            P_x_dists = [self.known_allele_contribs[pair]
                                                for pair in ind_genotype_pairs]

            return scipy.sparse.vstack(P_x_dists)

        # print "Calculating", len(unknown_pairs), "new dists"

        ## Get allele contribution hists and parent genotypes for unknown inds
        unknown_inds, unknown_genotypes = zip(*unknown_pairs)
        unknown_hist_indices = [self.ind_dict[ind] for ind in unknown_inds]
        unknown_contribs = self.all_contrib_hists[unknown_hist_indices]

        ## Normalize contributions of boundary inds, to convert counts to
        ## probabilities
        unknown_probs = normalize_sparse_array(unknown_contribs)

        ## Adjust contributions given that inds are in the boundary of the
        ## allele transmission tree, and may not receive an affected allele
        ## from their parents
        new_boundary_contribs = boundary_contribution(unknown_probs,
                                                            unknown_genotypes)

        ## Get the likelihood of each boundary ind having produced the
        ## observed allele frequency parameters, ie. P(k) given N, n
        boundary_control_liks = ind_control_liks(new_boundary_contribs,
                                                        self.N, self.n, self.k)

        ## Add result to known allele contribs
        for pair, row in zip(unknown_pairs, boundary_control_liks):
            self.known_allele_contribs[pair] = row

        ## Return P_x distribs for all inds
        P_x_dists = [self.known_allele_contribs[pair]
                                            for pair in ind_genotype_pairs]

        return scipy.sparse.vstack(P_x_dists)


def normalize_sparse_array(sparse_array):
    """
    Returns a normalized csr_matrix of the provided sparse array
    """
    coo_contribs = sparse_array.tocoo()
    row_sums = sparse_array.sum(axis=1)
    norm_data = np.array([d / row_sums[i]
                for i, d in izip(coo_contribs.row, coo_contribs.data)]).ravel()
    norm_contribs = scipy.sparse.csr_matrix((norm_data,
                                    (coo_contribs.row, coo_contribs.col)))

    return norm_contribs


def boundary_contribution(ind_contribs, parent_genotypes):
    """
    Returns the allele frequency distribution of each ind in ind_contribs
    that corresponds with their parents having the genotypes specified in
    `parent_genotypes`
    """
    ## DONE: Convert to lil_matrix before modifying entries +p1 id:166
    ## Convert to lil_matrix for more efficient manipulation of values
    boundary_contribs = ind_contribs.tolil().astype(np.float64)

    for i, g in enumerate(parent_genotypes):
        ## If tree ind is a het, their offspring have only a 50 percent
        ## chance of receiving an allele. If tree ind is a hom, we don't need
        ## to adjust the distribution, which was simulated with all hets
        ## NOTE: Underestimates for inds with both parents in tree +p3 id:118
        if np.max(g) == 1:
            boundary_contribs[i, 0] += 1.
            boundary_contribs[i] /= boundary_contribs[i].sum()

    return boundary_contribs.tocsr()


def ind_control_liks(ind_contribs, N, n, k):
    """
    Calculates the likelihood of observing `x` alleles, where x = range(0, k+1),
    in `n` samples from a population of size `N`, for each allele frequency
    distribution in sparse array `ind_contribs`.
    """
    ## TODO: Write tests for this function +p4 id:165
    ## Containers for new sparse array data
    newrow, newcol, newdata = [], [], []

    ## Calculate control liks for each row in ind_contribs
    for i, row in enumerate(ind_contribs):
        ## Get range of nonzero columns of each row, so we know how many
        ## values to iterate through
        min_col = np.min(row.tocoo().col)
        max_col = np.max(row.tocoo().col)
        alleles_range = range(min_col, np.min([max_col+1, N]))

        ## Find the probability of having drawn any number of alleles, up to
        ## the observed number, from the ind corresponding to this row
        for x in range(0, k+1):
            ## Iterate through possible real (unobserved) numbers of alleles in
            ## the whole population of size `N`, from this ind
            ##TODO: Case k == 0, sum row[0, K] only (weights are all 1) +t1
            P_x = np.sum([scipy.stats.hypergeom.pmf(x, N, K, n) * row[0, K]
                                                    for K in alleles_range])

            ## Append to new sparse array
            newrow.append(i)
            newcol.append(x)
            newdata.append(P_x)

    return scipy.sparse.csr_matrix((newdata, (newrow, newcol)))


def sum_control_liks(control_liks, k=None):
    """
    Returns the convolution of allele frequency distributions in control_liks.
    If the number of observed alleles `k` is provided, the result of each
    pairwise convolution is trimmed to `k` values, since larger values are
    inconsistent with the observation.
    """
    to_convolve = [row.toarray().ravel() for row in control_liks]

    while len(to_convolve) > 1:
        new_distribs = []

        ## Iterate through pairs of distributions
        for i in range(0, len(to_convolve)-1, 2):
            new_distrib = np.convolve(to_convolve[i], to_convolve[i+1])

            ## Truncate convolved distribution if `k` is specified
            if k is not None:
                new_distrib = new_distrib[:k+1]

            new_distribs.append(new_distrib)

        ## Update list of distribs to convolve
        to_convolve = new_distribs

    total = to_convolve[0]

    return total
