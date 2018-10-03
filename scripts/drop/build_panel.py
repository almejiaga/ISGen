# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 13:08:27 2015

@author: dominic
"""
from __future__ import division
import csv
from collections import defaultdict, Counter
import numpy as np
import copy
import random
import time
import os

class panel:
    def __init__(self, sim_results_path,proband_file, count_self = False):
        self.sim_results_path = sim_results_path
        self.clause_list = []
        assert count_self is True or count_self is False
        self.count_self = count_self
        p1 = time.time()

        with open(self.sim_results_path, 'rb') as sim_file:
            reader = csv.reader(sim_file,  delimiter=',' )
            self.probands = next(reader)
            self.Iterations = sum(1 for row in reader) / 2
            print str(self.Iterations) + ' iterations'
            indices = range(len(self.probands))
            self.proband_dict = dict(zip(self.probands, indices))
            ## Build dictionary of individuals who share the same allele,
            ## i.e. these are the clauses of individuals who cover each other
        p2 = time.time()

	    ## If you have a list of probands which is different from the list of
        ## all probands in the pedigree
        self.selected_probands = self.probands
        self.selected_probands_dict = dict(zip(self.selected_probands, indices))
        if proband_file != "None":
            self.selected_probands = map(str,[line.strip() for line in open(proband_file, 'rb')])
            indices = range(len(self.selected_probands))
            self.selected_probands_dict = dict(zip(self.selected_probands, indices))

        print "Finding allele clusters..."
        with open(self.sim_results_path, 'rb') as sim_file:
            header = sim_file.readline()
            while True:
                line1 = sim_file.readline()
                alleles1 = line1.split(',')
                line2 = sim_file.readline()
                alleles2 = line2.split(',')
                if not line2: break
            # for line in csv.reader(sim_file,  delimiter=','):
                build_allele_dict = defaultdict(list)
                for i, allele in enumerate(alleles1):
                    build_allele_dict[allele].append(self.probands[i])
                for i, allele in enumerate(alleles2):
                    build_allele_dict[allele].append(self.probands[i])
                ## Every clause from each simulation gets added to the list
                ## of all clauses. Alleles themselves are not used again.
                if self.count_self is False:
                    for value in build_allele_dict.values():
                        if len(value) > 1:
                            self.clause_list.append(value)
                elif self.count_self is True:
                    for value in build_allele_dict.values():
                        self.clause_list.append(value)

        p3 = time.time()
        print "Calculating proband coverage scores..."
        ## A list of coverages, in the same order as the list of probands.
        self.coverage = np.zeros(len(self.probands))
        self.coverage = list(self.coverage)
        self.selected_coverage = np.zeros(len(self.selected_probands))
        self.selected_coverage = list(self.selected_coverage)
        ## Ignore first dictionary since it contains the header
        for clause in self.clause_list:
            if self.count_self is True:
                ## Iterate through set of clause to avoid counting the clause
                ## twice for any homozygotes who may have both alleles in the
                ## clause.
                for ind in set(clause):
                    self.coverage[self.proband_dict[ind]] += len(clause)
                    try:
                        self.selected_coverage[self.selected_probands_dict[ind]] += len(clause)
                    except KeyError:
                        pass
            elif self.count_self is False:
                for ind in set(clause):
                    self.coverage[self.proband_dict[ind]] += len(clause) - 1
                    try:
                        self.selected_coverage[self.selected_probands_dict[ind]] += len(clause) - 1
                    except KeyError:
                        pass
            else:
                print "Invalid option! count_self must be True or False"
                sys.exit()

        p6 = time.time()
        ## A dictionary of individuals paired with the clauses which
        ## contain them.
        ##
        ## Each individual is only in a clause once, although if they have
        ## both alleles in a clause (ie they are homozygous), then they will
        ## both be included in the coverage score.
        ##
        ##@@ This could be quite slow for the number of clauses we will have
        ## in a large pedigree.
        print "Pairing probands and allele clusters..."
        self.ind_clause_dict = defaultdict(set)
        for i in range(len(self.clause_list)):
            for ind in self.clause_list[i]:
                self.ind_clause_dict[ind].add(i)

        p7 = time.time()

        self.panelbuildtimes = [p2 - p1, p3 - p2, p6 - p3, p7 - p6]


    def panel_probands_converage(self, panel):
        """
        Takes a list of individuals in panel, and outputs the probands they
        cover, along with the number of times they are covered, according to
        the allele dropping simulation file provided.
        """
        covered_probands = []
        clause_lengths = []
        for ind in panel:
            for clausenum in self.ind_clause_dict[ind]:
                clause_lengths.append(len(self.clause_list[clausenum]))
                covered_probands.extend(self.clause_list[clausenum])

        return Counter(covered_probands), clause_lengths


    def greedypanel(self, numinds):
        bestinds = []
        bestinds_coverage = []
        self.update_selected_coverage = copy.deepcopy(self.selected_coverage)
        ## Indicates whether clauses have been counted (0) or not (1)
        self.clause_status = np.ones(len(self.clause_list))
        # sel_ind_count = 0
        # count = 0
        #for run in range(numinds)
        runs = np.min([numinds, len(self.selected_probands)])
        for run in range(runs):
            ##@@ It would be fastest to save the coverage here rather than
            ## calculate it again.
            # bestind = self.probands[self.update_coverage.index(max(self.update_coverage))]
            bestind = self.selected_probands[self.update_selected_coverage.index(max(self.update_selected_coverage))]
            bestinds.append(bestind)
            bestinds_coverage.append(max(self.update_selected_coverage))

            for clausenum in self.ind_clause_dict[bestind]:
                ## Iterate through set of clause to avoid counting the clause
                ## twice for any homozygotes who may have both alleles in the
                ## clause.
                for ind in set(self.clause_list[clausenum]):
                    try:
                        if self.count_self is True:
                            self.update_selected_coverage[self.selected_probands_dict[ind]] -= self.clause_status[clausenum] * \
                                                                    (len(self.clause_list[clausenum]))
                        elif self.count_self is False:
                            self.update_selected_coverage[self.selected_probands_dict[ind]] -= self.clause_status[clausenum] * \
                                                                    (len(self.clause_list[clausenum]) - 1)
                        else:
                            print "Invalid option! count_self must be True or False"
                    except KeyError:
                        pass
                self.clause_status[clausenum] = 0

        return bestinds, bestinds_coverage


    def randompanel(self, numinds):
        ## Chooses random sample of indices, the individuals from those indices
        rand_indices = random.sample(xrange(len(self.selected_probands)), numinds)
        rand_inds = [self.selected_probands[i] for i in rand_indices]

        return rand_inds


    def getcoverage(self, indlist, proportion = False):
        self.clause_status = np.ones(len(self.clause_list))
        self.update_coverage = copy.deepcopy(self.coverage)
        panel_coverage = 0

        ind_coverage = []
        for panel_ind in indlist:
            panel_coverage += self.update_coverage[self.proband_dict[panel_ind]]
            ind_coverage.append(self.update_coverage[self.proband_dict[panel_ind]])

            for clausenum in self.ind_clause_dict[panel_ind]:
                ## Iterate through set of clause to avoid counting the clause
                ## twice for any homozygotes who may have both alleles in the
                ## clause.
                for ind in set(self.clause_list[clausenum]):
                    if self.count_self is True:
                        self.update_coverage[self.proband_dict[ind]] -= self.clause_status[clausenum] * \
                                                                    (len(self.clause_list[clausenum]))
                    elif self.count_self is False:
                        self.update_coverage[self.proband_dict[ind]] -= self.clause_status[clausenum] * \
                                                                    (len(self.clause_list[clausenum]) - 1)
                    else:
                        print "Invalid option! count_self must be True or False"
                        sys.exit()
                self.clause_status[clausenum] = 0

        if proportion is True:
            num_alleles = 2. * self.Iterations * len(self.probands)
            panel_coverage = 1. * panel_coverage / num_alleles

        return panel_coverage, ind_coverage


    def writefile(self, outfile, output):
        with open(outfile, 'a') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(output)


    def createfile(self, outfile):
        file_suffix = 1
        if os.path.isfile(outfile) is True:
            file, ext = os.path.splitext(outfile)
            while os.path.isfile(outfile) is True:
                file_suffix += 1
                outfile = file + '_' + str(file_suffix) + ext

            print "Warning: output file already exists. \
                                    Outputting to new file", outfile

        return outfile


    def run(self, panel_out_path, panelsizes, random_panels):
        sim_results_name, datafile_ext = os.path.splitext(os.path.basename(self.sim_results_path))
        outfile = os.path.join(panel_out_path, sim_results_name + "_panel.csv")

        newoutfile = self.createfile(outfile)

        runtimes = []
        greedypanels = []
        greedycoverage = []
        randomcoverage = []
        ind_coverage = []

        if panelsizes is None:
            maxinds = len(self.selected_probands)
            panelsizes = range(0, maxinds + 10, 10)

        for numinds in panelsizes:
            print "Building panels of", numinds, "individuals..."
            ## Build a panel of 'numinds' individuals to cover the pedigree
            a = time.time()
            greedypanel = self.greedypanel(numinds)
            greedypanels.append(greedypanel[0])

            ## Check coverage of panel
            b = time.time()
            cov, ind_cov = self.getcoverage(greedypanel[0])
            greedycoverage.append(cov)
            ind_coverage.append(ind_cov)

            # greedycoverage.append(self.getcoverage(greedypanel[0]))
            c = time.time()
            ## Check panel against random panel
            if random_panels != 0:
                randomtotal = []
                for i in range(random_panels):
                    d = time.time()
                    randompanel = self.randompanel(numinds)
                    e = time.time()
                    cov, ind_cov = self.getcoverage(randompanel)
                    randomtotal.append(cov)
                    # randomtotal.append(self.getcoverage(randompanel))
                    f = time.time()

                randomcoverage.append(np.mean(randomtotal))

                runtimes.append([b-a, c-b, e-d, f-e])
                self.writefile(newoutfile, [greedypanels[-1], greedycoverage[-1], randomcoverage[-1],
                                 runtimes[-1]])
                self.writefile(newoutfile, greedypanel[1])
            else:
                runtimes.append([b-a, c-b])
                self.writefile(newoutfile, [greedypanels[-1], greedycoverage[-1], 0,
                                     runtimes[-1]])
                self.writefile(newoutfile, greedypanel[1])

        return greedycoverage, randomcoverage, runtimes
