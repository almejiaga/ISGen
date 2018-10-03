from __future__ import division
import numpy as np
import sys, os
from collections import defaultdict, Counter, OrderedDict


class Pedigree:
    def __init__(self, pedfile):
        """
        Methods for analyzing pedigree files, with the columns specified
        in the header as:

            Ind    Father    Mother

        which are the only columns loaded, regardless of others present.
        """
        self.pedfile = pedfile
        alldata = np.genfromtxt(pedfile, skip_header=1, usecols=(0, 1, 2),
                                dtype=int)

        self.inds = alldata[:, 0]
        self.fathers = alldata[:, 1]
        self.mothers = alldata[:, 2]
        indices = range(len(alldata[:,0]))
        self.ind_dict = dict(zip(alldata[:,0], indices))

        ## Build a dictionary containing the parents and offspring of
        ## each individual
        self.parent_dict = dict(zip(self.inds, alldata[:, 1:3]))

        self.offspring_dict = defaultdict(list)
        for ind, father, mother in alldata:
            self.offspring_dict[father].append(ind)
            self.offspring_dict[mother].append(ind)

        ## Probands are individuals who are neither mothers nor fathers
        probands = set(self.inds).difference(set(self.fathers))
        self.probands = probands.difference(set(self.mothers))
        

    def ordered_lineage(self, ind):
        """
        Returns all ancestors of the given ind, along with the lengths of
        each lineage connecting to them
        """
        ordered_lineage = defaultdict(list)
        indnum = self.ind_dict[ind]
        current_generation = [self.mothers[indnum], self.fathers[indnum]]

        if current_generation[0] == current_generation[1] == 0:
            ordered_lineage[0].append(0)
            return ordered_lineage

        i = 0
        while len(current_generation) > 0:
            i += 1
            next_generation = []
            for ind in current_generation:
                ## If there are multiple paths, we save both lengths in the list
                ordered_lineage[ind].append(i)

            for ancestor in current_generation:
                indnum = self.ind_dict[ancestor]
                if self.mothers[indnum] != 0:
                    next_generation.append(self.mothers[indnum])
                if self.fathers[indnum] != 0:
                    next_generation.append(self.fathers[indnum])

            current_generation = next_generation

        return ordered_lineage


    def getlineage(self, ind):
        """
        Returns all the ancestors of an individual, including themselves,
        in a single list
        """
        lineage = [ind]
        father, mother = self.parent_dict[ind]

        if father != 0:
            lineage.extend(self.getlineage(father))

        if mother != 0:
            lineage.extend(self.getlineage(mother))

        return lineage
    

    def ordered_descendants(self, ind):
        """
        Returns a dictionary of all descendants of 'ind', with a corresponding
        list of path lengths from 'ind' to each descendant, in number of
        generations, as:

            descendants[ind] = [path_length_1, path_length_2, ... ]
        """
        descendants = defaultdict(list)
        descendants[ind] = [0]
        currentgen = self.offspring_dict[ind]

        i = 0
        while len(currentgen) > 0:
            i += 1
            nextgen = []
            for ind in currentgen:
                ## If there are multiple paths, we save all lengths in the list
                descendants[ind].append(i)
            for offspring in currentgen:
                nextgen.extend(self.offspring_dict[offspring])

            currentgen = nextgen

        return dict(descendants)


    def useful_ancestors(self, indlist):
        """Returns from the list of probands the list of ancestors shared
        by at least 2 probands. We keep only the points where coalescence could ever happend."""
        useful_ancs=[]
        
        for ind in indlist:
            for ind2 in indlist:
                if ind2 != ind:
                    common2 = []
                    common2.append(set(self.ordered_lineage(ind).keys()))
                    common2.append(set(self.ordered_lineage(ind2).keys()))
                    intersect2 = list(set.intersection(*common2))
                    
                    for a in intersect2:
                        if a not in useful_ancs:
                            useful_ancs.append(a)
            
        return useful_ancs


    def ancestors_dict(self, indlist):
        """Returns a dictionary lineages containing everyone in the probands
        ancestors. As values, in the lineages, we keep only the 'useful people' from
        the precedent function. Meaning that we store in the lineages only the
        possible coalescence point between at least 2 people."""
        lineages = dict()
        useful_ancs = set(self.useful_ancestors(indlist))
        
        for ind in indlist:
            ancs = self.ordered_ancestors(ind)
            ## ancestors distance when the distance is < N=10 
            used_ancs = useful_ancs.intersection(set(ancs.keys()))
            
            lin = dict()
            for a in used_ancs:
                lin[a] = ancs[a]
                
            lineages[ind] = lin
            
            for a2 in self.ordered_lineage(ind):
                ## We need to have every person that we could climb to
                ## in the dictionary. This is why we add everyone in 
                ## the lineage of ind.
                ancs_a = self.ordered_ancestors(a2)
                used_ancs_a = useful_ancs.intersection(set(ancs_a.keys()))
                lin_a = dict()
                for aa in used_ancs_a:
                    lin_a[aa] = ancs_a[aa]
                    
                lineages[a2] = lin_a
    
        return lineages

    
    def ordered_ancestors(self, ind, N=10):
        """
        Returns the list [ancestor, generation] for each ancestor of one individual.
        If there are multiple paths to reach one ancestor, the generation number bewteen the 2 individuals
        corresponds to the shortest path. Plus, we keep only the people that are closer
        than N=10 (by default) generations.
        """
        temp_ancs = dict()
        lineage = self.ordered_lineage(ind)

        for a in lineage:
            a_gen = min(lineage[a])
            if a_gen < N:
                temp_ancs[a] = a_gen

        return temp_ancs


    def allowedinds(self, indlist):
        """
        Returns:

            commonancestors - list of all common ancestors between inds in
                                indlist
            cone_inds - dict with format {anc: set(descendants of anc)}
            ind_cones - dict with format {ind: set(ancs ancestral to ind)}
        """
        ancestorsets = []
        for ind in indlist:
            ancestorsets.append(set(self.getlineage(ind)))

        commonancestors = list(set.intersection(*ancestorsets))

        if len(commonancestors) == 0:
            print "Warning: No common ancestors!"

        cone_inds = {}
        ind_cones = defaultdict(set)
        for anc in commonancestors:
            tmp_dict = self.ordered_descendants(anc)
            tmp = set(tmp_dict.keys())
            ## We only need to track the individuals who could possibly climb
            ## to some point on the tree. To find these individuals, we take
            ## the descendants of one common ancestor, and take the union of
            ## its intersection with the lineage of every affected proband.
            ## repeat for every common ancestor to get the descent cones.
            tmp_intersect = set()
            for ancestorset in ancestorsets:
                tmp_intersect.update(tmp.intersection(ancestorset))
            tmp = tmp_intersect
            cone_inds[anc] = tmp

            ## Build dict of cone memberships for each individual
            for ind in tmp:
                ind_cones[ind].add(anc)

            ## Nonexistant inds, denoted by '0', have no descendats
            ind_cones[0] = set()

        return commonancestors, cone_inds, ind_cones
