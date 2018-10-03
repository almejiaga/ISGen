from __future__ import division
import numpy as np
import random
import sys, os
from collections import defaultdict
from profilehooks import profile


class VPrint(object):
    def __init__(self):
        """ Easily turn print output on/off """
        self.verbose = False

    def __call__(self, *args):
        if self.verbose:
            # Print each argument separately so caller doesn't need to
            # stuff everything to be printed into a single string
            for arg in args:
                print (arg)

            print('')

## Instantiate printing class globally
vprint = VPrint()


class AncFinder:
    def __init__(self, ind_cones, parent_dict, lineages, args):
        """
        Methods for inferring allele transmission histories for alleles
        which share a single origin within a known genealogy
        """
        self.ind_cones = ind_cones
        self.parent_dict = parent_dict
        self.lineages = lineages
        self.args = args
        self.inheritance_dict = defaultdict(list)
        self.signals = dict()
        self.common_ancs = []
        self.genotypes = dict()
        
    def check_ind(self, ind):
        """
        Convenience function for checking coalescence possibilities from 'ind'
        """
        return len(self.ind_cones[ind].intersection(self.common_ancs))


    def is_founder(self, ind):
        """ Checks if ind is a founder, ie. both parents given as '0'. """
        return all(self.parent_dict[ind] == [0, 0])


    def hom_ancs(self, ind):
        """
        Check if the given individual has a path to coalescence through both
        parents, in which case they can be simulated as a homozygote.
        """
        ## TODO: Should we always coalesce in founders? +p1  id:178
        ## Founders are a special case since they have no parents. Their
        ## ability to be homozygotes does not depend on ancestors other
        ## than themselves.
        if self.is_founder(ind):
            return set([ind])

        gf, gm = self.parent_dict[ind]
        gp_common_ancs = self.ind_cones[gf].intersection(self.ind_cones[gm])

        return gp_common_ancs.intersection(self.common_ancs)


    def find_common_ancs(self, cur_carriers):
        """ Find the common ancestors among current carriers """
        all_ancs = [self.ind_cones[ind] for ind in cur_carriers]
        common_ancs = set.intersection(*all_ancs)

        return common_ancs


    def initial_signals(self,ind):
        """Initialisation of the signals"""
        """Sends the first signals to the probands ancestors"""
        self.signals[ind] += 1

        for ind_a in self.lineages[ind]:
            if ind_a != 0:
                self.signals[ind_a] += pow(0.5, int(self.lineages[ind][ind_a]))


    def update_signals(self,ind,parent,founders):
        """
        Updating the signals,
        step 1 : erasing the signals sent by ind
        step 2 : parent sends its own signals

        NOTE: if parent has already received an allele, he has already sent
        his signal
        """
        #self.signals[ind]-=1
        
        if self.signals[ind] < 0:
            self.signals[ind] = 0

        for ind_a in self.lineages[ind]:
            if ind_a != 0:
                self.signals[ind_a] -= pow(0.5, int(self.lineages[ind][ind_a]))

                if self.signals[ind_a] < 0:
                    self.signals[ind_a] = 0

        if self.genotypes[parent] == 0:
            self.signals[parent] += 1

            if parent not in founders:
                for parent_a in self.lineages[parent]:
                    if parent_a != 0:
                        self.signals[parent_a]+= pow(0.5,
                                self.lineages[parent][parent_a])


    def sum_parent_signals_cone(self, ind, parent):
        """
        Sum of the signals in the parent's cone without the signal sent by
        ind, The signals are balanced by the distance between the ancestor and
        the climbing individual.
        """
        s = 0
        s += self.signals[parent] - pow(0.5, 1) # TODO: Check why this is done
        for parent_a in self.lineages[parent]:
            if parent_a != 0:
                parent_gen = self.lineages[parent][parent_a]
                if parent_a in self.lineages[ind]:
                    ind_gen = self.lineages[ind][parent_a]
                    s += pow(0.5, parent_gen) * \
                            (self.signals[parent_a] - \
                            pow(0.5, ind_gen))

        if s < 0:
            vprint("Negative signal!")
            s = 0

        return s


    def add_genotypes(self, ind):
        """
        Updates the genotype of 'ind', who has been climbed to. If ind is not a
        carrier, they receive an allele. If ind is a heterozygous carrier, the
        new allele either coalesces, or ind becomes a homozygote if both
        alleles can still coalesce. If ind is a homozygote, nothing happens.
        """
        parent_alleles = self.genotypes[ind]

        ## Most climbing scenarios (counted by number of possible
        ## configurations, not most commonly encountered) do not involve
        ## continuing to climb through the parent. Set this as the
        ## default behaviour
        to_climb = None

        ## Initialize importance sampling likelihood
        lik = 1.

        ## If chosen parent is not a carrier yet, they receive an allele
        ## and continue climbing
        if parent_alleles == 0:
            self.genotypes[ind] += 1
            to_climb = ind

        ## If chosen parent is a carrier, check if homozygote can be created
        elif parent_alleles == 1:
            ## If the existing allele is the only valid path, then we
            ## use importance sampling to force coalescence, and note that
            ## we don't have to continue climbing this path
            gp_common_ancs = self.hom_ancs(ind)
            if len(gp_common_ancs) == 0:
                lik *= 0.5

            ## Otherwise we have a chance to create a valid homozygote, in
            ## which case we continue to climb the new allele. Otherwise
            ## we're done with this part of the lineage and do nothing
            else:
                if np.random.random() < self.args.sim_homs:
                    ## Adjust for importance sampling favouring/not favouring
                    ## the creation of homozygotes
                    lik *= 0.5 / self.args.sim_homs
                    to_climb = ind
                    self.genotypes[ind] += 1

                    ## Adjust common ancestors to include those of both
                    ## grandparents, who must have both been carriers
                    self.common_ancs = self.common_ancs.intersection(gp_common_ancs)

                else:
                    ## Adjust for importance sampling favouring/not favouring
                    ## the creation of homozygotes
                    lik *= 0.5 / (1 - self.args.sim_homs)

        ## If chosen parent is a homozygote, we're done - no need to check
        return lik, to_climb
    
    
    def climb_early(self, ind):
        """
        Chooses a parent allele to climb to from the given allele, ensuring
        that the choice is consistent with coalescence with the other alleles.
        Here, we try to chose the parent that is the most likely to give early
        coalescence.
        """
        ## Initialize importance sampling likelihood
        lik = 1.

        ## Check which parents can be climbed to and still allow all alleles
        ## to coalesce, ensuring only one allele is inherited from each parent
        parents = self.parent_dict[ind]
        can_inherit = [p for p in parents if p not in self.inheritance_dict[ind]]
        shared_ancs = [self.check_ind(p) for p in can_inherit]

        valid_parents = []
        num_anc_weights = []
        signal_weights = []
        for c, s in zip(can_inherit, shared_ancs):
            if s > 0:
                valid_parents.append(c)
                num_anc_weights.append(s)
                if self.args.kinship:
                    signal_weights.append(self.sum_parent_signals_cone(ind,c))

        if self.args.kinship:
            sample_weights = signal_weights
        else:
            sample_weights = num_anc_weights

        ## Calculate importance sampling factor based on number of valid
        ## parents. This will only introduce such a factor if one of two
        ## parents is the sole direction consistent with coalescence.
        if len(can_inherit) == 2 and len(valid_parents) == 1:
            lik *= 0.5
            climb_parent = valid_parents[0]
            
        elif len(can_inherit) == 2 and len(valid_parents) == 2:
            # for c in valid_parents:
            #     signal_weights.append(self.sum_parent_signals_cone(ind,c))
            if np.sum(sample_weights) != 0:
                sample_weights = np.array(sample_weights) / np.sum(sample_weights)
            else:
                sample_weights=[0.5, 0.5]
                
            ## Weight parent choice by number of shared ancs with other
            ## carriers
            climb_index = np.random.choice(range(len(valid_parents)),
                                            p=sample_weights)
            climb_parent = valid_parents[climb_index]
            lik *= 0.5 / sample_weights[climb_index]
            self.inheritance_dict[ind].append(climb_parent)
            
        elif len(can_inherit) == 1 and len(valid_parents) == 1:
            climb_parent = valid_parents[0]
        elif len(valid_parents) == 0:
            assert self.is_founder(ind)
            climb_parent = ind

        return lik, climb_parent
    

    def coalesce(self, genotypes, signals):
        """
        Climbs known homozygotes and heterozygotes until they coalesce in
        a common ancestor.
        """
        ## Initialize likelihood
        loglik = 0.

        ## Track inheritance of each allele
        self.signals = signals
        self.genotypes = genotypes
        self.inheritance_dict = defaultdict(list)

        ## Find inds to climb and common ancestors
        inds_to_climb = genotypes.keys()
        founders = set()
        self.common_ancs = self.find_common_ancs(inds_to_climb)

        for ind in inds_to_climb:
            self.initial_signals(ind)
        
        ## Climb alleles until they coalesce
        while len(inds_to_climb) > 0:
            ## Individuals to climb next
            climb_next = []

            ## Go through inds to climb in random order
            for ind in random.sample(inds_to_climb, len(inds_to_climb)):
                ## Climb each allele in the current individual
                steplik, parent = self.climb_early(ind)
                self.update_signals(ind, parent,founders)
                
                ## Update genotyped of the individual climbed to
                innerlik, newclimber = self.add_genotypes(parent)
                
                ## Update overall likelihood from this climbing step
                loglik += np.log2(steplik) + np.log2(innerlik)

                ## Update list of common ancestors and carriers to climb,
                ## noting that founders are denoted by having '0' for both
                ## parents
                if newclimber is not None:
                    new_ancs = self.ind_cones[newclimber]
                    self.common_ancs = self.common_ancs.intersection(new_ancs)

                    ## Check if we've climbed to a founder
                    if self.is_founder(newclimber):
                        founders.add(newclimber)
                    else:
                        climb_next.append(newclimber)

            inds_to_climb = climb_next

        ## Make sure we've only selected one founder
        assert len(founders) == 1

        return next(iter(founders)), loglik, genotypes
