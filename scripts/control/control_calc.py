# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 14:10:43 2015

@author: dominic
"""
from __future__ import division
from functools import reduce
import numpy as np
import scipy.stats
import scipy.sparse
import sys
import csv
from collections import defaultdict, OrderedDict, Counter
import tables
from profilehooks import timecall


def writefile(output, outfile, write_mode):
    with open(outfile, write_mode) as outfile:
        writer = csv.writer(outfile)
        writer.writerow(output)


def tree_mean_var(tree, mean_var_dict):
    """
    Returns the mean and variance of the contribution of the tree, using
    the dict mean_var_dict[ind] = [ind_mean, ind_variance]

    The specific patients used to generate the tree are excluded from the
    contribution.
    """
    tot_mean = 0
    tot_squared_mean = 0
    tot_var = 0
    for ind in self.getboundary(tree):
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


class Pedigree:
    def __init__(self, pedfile):
        alldata = np.genfromtxt(fname = pedfile,
                                skip_header = 1,
                                usecols = (0, 1, 2))

        indices = range(len(alldata[:,0]))
        self.inds = map(int, alldata[:,0].tolist())
        fatherlist = map(int, alldata[:,1].tolist())
        self.motherlist = map(int, alldata[:,2].tolist())
        self.ind_dict = dict(zip(alldata[:,0], indices))
        self.father_dict = dict(zip(map(int,alldata[:,0]),map(int,alldata[:,1])))
        self.mother_dict = dict(zip(map(int,alldata[:,0]),map(int,alldata[:,2])))
        self.offspring_dict = defaultdict(list)

        for i, ind in enumerate(self.motherlist):
            self.offspring_dict[ind].append(self.inds[i])
        for i, ind in enumerate(fatherlist):
            self.offspring_dict[ind].append(self.inds[i])

        self.probands = self.getprobands()


    def getprobands(self):
        ## Probands are individuals who are not mothers or fathers
        motherset = set(self.mother_dict.values())
        fatherset = set(self.father_dict.values())
        probands = [self.inds[i] for i,n in enumerate(self.inds) if
                    (self.inds[i] not in motherset and
                    self.inds[i] not in fatherset)]

        return set(probands)


    def getleaves(self):
        ## Leaves are those individuals without parents, indicated by a '0' in
        ## the pedigree.
        leaves=[self.inds[i] for i,n in enumerate(self.motherlist) if n==0]
        return leaves


    def ordered_descendants(self, ind, out_type = 'dict'):
        ordered_descendants = defaultdict(list)
        ordered_descendants[ind] = [0]
        current_generation = self.offspring_dict[ind]

        i = 0
        while len(current_generation) > 0:
            i += 1
            next_generation = []
            for ind in current_generation:
                ## If there are multiple paths, we add both lengths in the list
                ordered_descendants[ind].append(i)
            for offspring in current_generation:
                next_generation.extend(self.offspring_dict[offspring])

            current_generation = next_generation

        if out_type == 'dict':
            return dict(ordered_descendants)
        elif out_type == 'set':
            return set(ordered_descendants.keys())
        elif out_type == 'list':
            return ordered_descendants.keys()


    def buildpath(self, startlist=None):
        """
        Builds a path through the pedigree, in the form of lists of individuals
        who can receive an allele from at least one parent at each generation.

        Note this is different from a typical allele dropper, where inds
        receive an allele from each parent, and can have only two. We are
        interested in the shortest path an allele could possibly take to each
        proband, so we drop all of them with probability 1.

        Thus the path is built based on a single parent having an allele, not
        both as would normally be expected. Each allele is passed immediately,
        then deleted from the parent.
        """
        ## Unless provided with a starting list, we start with the founders
        ## (leaves) of the pedigree
        ##@@ This is broken... shouldn't be an option
        if startlist is None:
            startlist = self.getleaves()

        hasalleles = defaultdict(int) # 0:no alleles, 1:has alleles
        path = [] # To store the individuals to receive alleles at each step
        nextgen= [] # Stores the individuals who can receive alleles next
        nextgen.extend(startlist) # First step is the founders

        ## Individuals for next step are chosen from offspring of previous step
        ## Offspring with only one allele-calculated parent get carried over.
        offspring_notselected = []

        while len(nextgen) > 0:
            ## Add step to path
            path.append(nextgen)

            ## Assign alleles to current step
            for ind in nextgen:
                hasalleles[ind] = 1

            ## Find offspring of current step
            offspring = []
            offspring.extend(offspring_notselected)
            offspring_notselected = []
            for ind in nextgen:
                offspring.extend(self.offspring_dict[ind])

            ##@@ Could use a set instead? Do we need to preserve order?
            offspring = list(OrderedDict.fromkeys(offspring))

            ## Find next step from offspring of current step.
            nextgen = []
            for off in offspring:
                ##@@ I think some of these checks are unnecessary
                if (hasalleles[off] == 0 and
                        self.mother_dict[off] != 0 and
                        self.father_dict[off] != 0 and
                        (hasalleles[self.mother_dict[off]]==1 and
                        hasalleles[self.father_dict[off]]==1)):
                    ## We append individuals whose parents have alleles, and
                    ## remove their parents. This is done so the parents
                    ## can then inherit another allele later.
                    nextgen.append(off)
                else:
                    offspring_notselected.append(off)

        return path


    def passalleles(self, path, homhet = 'het'):
        """
        Drop alleles to simulate contribution to probands from every individual
        in the pedigree.

        We drop from the founders, but also track all alleles as they
        pass through descendants, so we can get statistics on all individuals
        in the pedigree at once.
        """
        numremaining = len(self.inds)
        ## Initialize alleles as empty lists. Start with leaves/ancestors,
        ## whose alleles we are tracking.
        alleles = defaultdict(list)
        startlist = self.getleaves()
        numremaining -= len(startlist)

        ## Alleles for the leaves/ancestors are simply their own integer label.
        for ind in startlist:
            if homhet == 'hom':
                alleles[ind] = [[str(ind) + 'a'], [str(ind) + 'b']]
            elif homhet == 'het':
                alleles[ind] = [[0], [str(ind) + 'b']]

        ## Each time we complete a loop, the alleles of more individuals can be
        ## computed. We continue until all individuals have alleles and the
        ## nextgen list is empty.
        for i in range(1, len(path)):
            numremaining -= len(path[i])
            for ind in path[i]:
                ## Mark if an individual is in the tree provided, if any.
                if homhet == 'hom':
                    alleles[ind] = [[str(ind) + 'a'], [str(ind) + 'b']]
                elif homhet == 'het':
                    alleles[ind] = [[0], [str(ind) + 'b']]
                mother = self.mother_dict[ind]
                father = self.father_dict[ind]

                ## Each individual gets one allele from each parent, each
                ## passed down with a 50% chance.
                passallele1 = alleles[mother][np.random.randint(2)]
                passallele2 = alleles[father][np.random.randint(2)]
                alleles[ind][0].extend(passallele1)
                alleles[ind][1].extend(passallele2)

        return alleles


    def get_contrib(self, allele_dict, region_dict=None):
        """
        Takes a dictionary of all alleles received by probands of the pedigree
        in one allele drop simulation, and returns the number of copies of
        each allele received.

        Individuals all pass one of their two alleles to each of their
        offspring.

        If a region dict is provided to sort probands by region, the function
        returns a dictionary of allele counts per individual per region.
        """
        if region_dict is None:
            region_dict = {'All Probands': self.probands}
        else:
            ## Any proband without a specified region gets assigned 'None'
            hasregion = reduce(set.union, map(set, region_dict.values()))
            noregion = self.probands.difference(hasregion)

            ## Don't count inds in automatically-generated 'All Probands'
            num_inds_given_regions = np.sum([len(region_dict[r])
                        for r in region_dict.keys() if r != 'All Probands'])

            if num_inds_given_regions > len(hasregion):
                print "------------------------------------------"
                print "Warning: individuals with multiple regions"
                print "------------------------------------------"

            if len(noregion) > 0:
                region_dict['No Region'] = noregion

        all_alleles = defaultdict(list)
        regions = sorted(region_dict.keys())

        ## Collect the alleles from the probands in each region
        for region, prob_list in region_dict.iteritems():
            for prob in prob_list:
                ## If we have sub-sampled the pedigree, not all probands
                ## will be present
                if prob in allele_dict:
                    all_alleles[region].extend(allele_dict[prob][0])
                    all_alleles[region].extend(allele_dict[prob][1])

        contribs = scipy.sparse.lil_matrix((len(self.inds),
                                                    len(region_dict.keys())))

        ## Add up contributions for each ancestor by region.
        for i in range(len(regions)):
            alleles = all_alleles[regions[i]]
            anc_counts = Counter([int(x[0:-1]) for x in alleles
                                                            if type(x) == str])

            for j in range(len(self.inds)):
                contribs[j, i] += anc_counts[self.inds[j]]

        return contribs


def load_contrib(contrib_file):
    alldata = np.genfromtxt(fname = contrib_file,
                            skip_header = 0,
                            delimiter = ',',
                            usecols = (0, 1, 2))

    mean_var_dict = {}
    for row in alldata:
        mean_var_dict[int(row[0])] = map(float, row[1:])

    return mean_var_dict


class Tree(Pedigree):
    def __init__(self, pedfile, treefile):
        """
        Functions to handle analysis of the allele inheritance trees contained
        in treefile. Inherits from Pedigee.
        """
        Pedigree.__init__(self, pedfile)
        self.treefile = treefile
        with tables.open_file(treefile, 'r') as f:
            self.alltrees = f.root.trees[:]
            self.allgenotypes = f.root.genotypes[:]
            self.climb_liks = map(lambda x: 2.**x, f.root.liks[:])


    def getboundary(self, tree_num):
        """
        Returns the individuals who are the offspring of an individual in the
        tree, but not in the tree themselves. Also returns the genotype of
        the parents of these individuals
        """
        ## Load trees and their associated genotypes
        tree = self.alltrees[tree_num]
        genotypes = self.allgenotypes[tree_num]

        ## Use tuples so that the genotypes are hashable
        boundary_inds = defaultdict(tuple)

        ## Iterate through tree inds and their genotypes
        for ind, genotype in zip(tree, genotypes):
            for offspring in self.offspring_dict[ind]:

                ## Don't include individuals in the tree
                if offspring not in set(tree):
                    ## Use temporary list to 'append' to tuple
                    genotype_list = list(boundary_inds[offspring])
                    genotype_list.append(genotype)
                    boundary_inds[offspring] = tuple(genotype_list)

        ## Make sure only two parents have been added per individual
        assert len(boundary_inds) > 0, ("No individuals in tree boundary!",
                                                        self.alltrees[tree_num])
        assert np.max([len(b) for b in boundary_inds.values()]) <= 2

        return boundary_inds


    def sort_tree_stats_by_anc(self, mean_var_dict):
        tree_stats_by_anc = defaultdict(list)
        for tree in self.alltrees:
            anc = self.find_founder(tree)
            mean,var = self.tree_mean_var(tree, mean_var_dict)
            tree_stats_by_anc[anc].append((mean,var))

        return tree_stats_by_anc


    def passalleles_tree(self, treelist, homhet = 'het'):
        """
        Drop affected alleles from the provided tree, to get the distribution
        in the probands.

        Returns a dict of all individuals in the pedigree and their alleles.
        """
        ## Initialize alleles as empty lists. Start with leaves/ancestors,
        ## whose alleles we are tracking.
        alleles = defaultdict(lambda: [0,0])
        nextgen = treelist
        tree = set(treelist)

        ## Alleles for the leaves/ancestors are simply their own integer label.
        for ind in treelist:
            if homhet == 'het':
                alleles[ind] = [0,1]
            elif homhet == 'hom':
                alleles[ind] = [1,1]
            else:
                print "Invalid homhet option!"
                sys.exit()

        ## Each time we complete a loop, the alleles of more individuals can be
        ## computed. We continue until all individuals have alleles and the
        ## nextgen list is empty.
        i = 0
        while len(nextgen) > 0:
            i += 1
            ## Find offspring of current step
            offspring = []
            for ind in nextgen:
                offspring.extend(self.offspring_dict[ind])

            offspring = list(set(offspring).difference(tree))

            for ind in offspring:
                mother = self.mother_dict[ind]
                father = self.father_dict[ind]
                passallele1 = alleles[mother][np.random.randint(2)]
                passallele2 = alleles[father][np.random.randint(2)]

                ## Pass alleles to offspring
                alleles[ind][0] = passallele1
                alleles[ind][1] = passallele2

            nextgen = offspring

        return alleles


    def sim_tree_freq(self, tree, homhet = 'het'):
        """
        Simulates the allele frequency on the given tree. Individuals are
        either homozygous or heterozygous carriers depending on how 'homhet'
        is set
        """
        alleles = self.passalleles_tree(tree, homhet)
        probs_notree = self.probands.difference(set(tree))

        allele_count = 0
        for prob in probs_notree:
            allele_count += np.sum(alleles[prob])

        return allele_count


    def freq_to_carriers(self, freq, numprobs):
        """
        Returns the number of carriers that would be expected at frequency freq
        """
        return int(np.round(freq * numprobs, 0))


    def find_founder(self, tree):
        """
        Returns the founder of the tree. Checks every individual even after a
        founder has been detected, to handle errors that may have produced a
        tree with multiple founders.
        """
        founders = []
        for ind in tree:
            if (self.mother_dict[ind] == 0 and
                    self.father_dict[ind] == 0):
                founders.append(ind)

        if len(founders) != 1:
            print "Founders:", founders
            raise ValueError("Error: tree should have exactly one founder.")

        return founders[0]
