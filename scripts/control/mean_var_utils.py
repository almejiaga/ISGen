import numpy as np


def tree_mean_var(boundary_inds, mean_var_dict):
    """
    Returns the mean and variance of the contribution of the tree, using
    the dict mean_var_dict[ind] = [ind_mean, ind_variance]

    The specific patients used to generate the tree are excluded from the
    contribution.
    """
    tot_mean = 0
    tot_squared_mean = 0
    tot_var = 0
    for ind in boundary_inds:
        ## The boundary does not include tree inds, so the mean and
        ## variance excludes the specific patients who generated the tree.
        mean, var = mean_var_dict[ind]
        tot_mean += mean
        tot_squared_mean += mean ** 2
        tot_var += var

    tree_mean = 0.5 * tot_mean
    tree_var = 0.5 * tot_var + 0.25 * tot_squared_mean

    return tree_mean, tree_var


def tree_correction(tree_mean, tree_var, freq_min, freq_max):
    if tree_var == 0:
        ## If the variance is zero, we return 1 if the mean is in the
        ## provided range, and zero otherwise.
        return 1 * (tree_mean >= freq_min and tree_mean <= freq_max)

    left = scipy.stats.norm(tree_mean, np.sqrt(tree_var)).cdf(freq_min)
    right = scipy.stats.norm(tree_mean, np.sqrt(tree_var)).cdf(freq_max)

    return right - left


