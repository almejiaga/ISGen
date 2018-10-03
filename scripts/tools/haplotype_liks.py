import tables
import numpy as np
import sys, os
import scipy.stats
import argparse

sys.path.append(os.path.expanduser('~/project/anc_finder/scripts/climb/'))
import ped


def get_num_hidden_transissions(Pedigree, tree):
    founders = []
    tree = set(tree)
    for ind in tree:
        ind_num = Pedigree.ind_dict[ind]
        if Pedigree.mothers[ind_num] == Pedigree.fathers[ind_num] == 0:
            founders.append(ind)
            
    assert len(founders) == 1
    ind = founders[0]
    offspring_in_tree = tree.intersection(Pedigree.offspring_dict[ind])
    
    ## The first hidden transmission was when the founder received the allele
    num_hidden_transmissions = 1
    
    while len(offspring_in_tree) == 1:
        num_hidden_transmissions += 1
        ind = next(iter(offspring_in_tree))
        offspring_in_tree = tree.intersection(Pedigree.offspring_dict[ind])
            
    return num_hidden_transmissions


def main(args):
    pedfile = os.path.expanduser(args.pedfile)
    climb_lik_file = os.path.expanduser(args.climb_lik_file)
    control_lik_file = os.path.expanduser(args.control_lik_file)
    hap_length = args.haplotype_length

    print "Loading pedigree..."
    P = ped.Pedigree(pedfile)

    with tables.open_file(climb_lik_file, 'r') as f:
        shape = f.root.liks.shape
        haplotype_liks = np.zeros(shape)

        for i, tree in enumerate(f.root.trees):

            num_hidden_transmissions = get_num_hidden_transissions(P, tree)
            scale = 1. / (len(tree) - num_hidden_transmissions)
            haplotype_liks[i] = scipy.stats.erlang(2, scale=scale).pdf(hap_length)

            if i % 10000 == 0:
                print "Calculating haplotype likelihoods for tree number", i, \
                      "of", len(haplotype_liks)
                
    norm_hap_liks = np.log2(haplotype_liks / np.sum(haplotype_liks))

    print "Writing combined likelihoods to file..."
    with tables.open_file(control_lik_file, 'a') as hapfile:
        tot_hap_liks = hapfile.create_carray(hapfile.root, 'tot_hap_liks',
                    tables.FloatAtom(), shape=shape,
                    title='Total tree likelihoods, including haplotype length')
        haps = hapfile.create_carray(hapfile.root, 'haplotype_liks',
                tables.FloatAtom(), shape=shape,
                title='Likelihood of observed haplotype length:' + \
                        str(hap_length) + 'Morgans')

        haps[:] = norm_hap_liks
        tot_liks = hapfile.root.tot_liks[:]
        t = norm_hap_liks + tot_liks
        tot_hap_liks[:] = norm_hap_liks + tot_liks


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("pedfile",
            help="Pedigree in which climbing simulations were performed")
    parser.add_argument("climb_lik_file",
            help="File containing output of climbing simulations")
    parser.add_argument("control_lik_file",
            help="File containing output of control likelihood calculations")
    parser.add_argument("haplotype_length", type=float,
            help="Length in Morgans of the shared haplotype observed in " +\
                    "the known patients and carriers")

    args = parser.parse_args()
    main(args)

