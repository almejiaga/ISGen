import numpy as np
import tables
import sys, os


def init_array(outfile, array_name):
    """ Initialize extendable arrays to store control likelihoods """
    with tables.open_file(outfile, 'w') as f:
        title = 'Likelihood of tree having produced the observed allele freq'
        f.create_earray(f.root, array_name, atom=tables.FloatAtom(),
                                title=title, shape=(0,))


def incremental_write(liks, trees, outfile):
    ## Placeholder for ancestor values
    ancs = np.zeros(len(liks))

    ## Write all output to hdf5 file, in separate nodes
    with tables.open_file(outfile, 'a') as f:
        f.root.liks.append(liks)
	f.root.ancs.append(ancs)

        ## Trees, which are variable-length, must be added individually
        for i, tree in enumerate(trees):
	    if i % 100 == 0:
		print "tree", i
	    tree = [int(t) for t in tree.split(',')]
            genotypes = np.ones(len(tree))
            f.root.trees.append(tree)
            f.root.genotypes.append(genotypes)


cor_lik_file = '/RQusagers/dnelson/project/anc_finder/results/BALSAC/CAID_3M_all_cor_liks.txt'
tree_anc_file = '/RQusagers/dnelson/project/anc_finder/results/BALSAC/CAID_3M_all_tree_anc.csv'
climb_lik_file = '/RQusagers/dnelson/project/anc_finder/results/BALSAC/CAID_3M_all_anc_out.csv'


max_trees = 1000

print "Loading climbing likelihoods"
climb_liks = np.genfromtxt(climb_lik_file, skip_header=True,
				delimiter=',', usecols=[1])[:max_trees]
print "Loading trees"
trees = [line.strip() for line in open(tree_anc_file, 'r')][:max_trees]

climb_outfile = os.path.expanduser('~/temp/CAID_climb_1000.h5')
with tables.open_file(climb_outfile, 'w') as f:
    ## Create extendable arrays so we can incrementally write output
    f.create_earray(f.root, 'liks', atom=tables.FloatAtom(), shape=(0,))
    f.create_earray(f.root, 'ancs', atom=tables.FloatAtom(), shape=(0,))

    ## Trees, which are variable-length, must be added individually
    f.create_vlarray(f.root, 'trees', atom=tables.IntAtom())
    f.create_vlarray(f.root, 'genotypes', atom=tables.IntAtom())

incremental_write(climb_liks, trees, climb_outfile)

## Store control likelihoods
# control_outfile = os.path.expanduser('~/temp/CAID_control.h5')
# 
# init_array(control_outfile, array_name='control_liks')
# with tables.open_file(control_outfile, 'a') as f:
#     f.root.control_liks.append([np.log2(conv_prob)])

