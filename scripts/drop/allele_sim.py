# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 10:29:23 2015

@author: Dom
"""

import numpy as np
from numpy.random import randint
import random
import os.path
import csv
from collections import defaultdict, OrderedDict

#import datetime
#import time
#import PyPedal


def loadpedfile(datafile):
    alldata = np.genfromtxt(fname = datafile,
                            skip_header = 1,
                            usecols = (0, 1, 2))

    indices = range(len(alldata[:,0]))
    inds = map(int, alldata[:,0].tolist())
    fathers = map(int, alldata[:,1].tolist())
    mothers = map(int, alldata[:,2].tolist())
    ind_dict = dict(zip(alldata[:,0], indices))
    father_dict = dict(zip(alldata[:,1], indices))
    mother_dict = dict(zip(alldata[:,2], indices))
    build_offspr_dict = defaultdict(list)

    for i, ind in enumerate(mothers):
        build_offspr_dict[ind].append(inds[i])
    for i, ind in enumerate(fathers):
        build_offspr_dict[ind].append(inds[i])
    offspr_dict = dict(build_offspr_dict)

    datalist = []
    datalist.extend([inds, fathers, mothers, ind_dict,
                     father_dict, mother_dict, offspr_dict])

    return datalist


def writefile(output, outfile, write_mode = 'wb'):
    with open(outfile, write_mode) as outfile:
        writer = csv.writer(outfile)
        writer.writerow(output)


def createfile(outfile, header_list, newfile = True):
    file_suffix = 1
    if os.path.isfile(outfile) is True and newfile is True:
        file, ext = os.path.splitext(outfile)
        while os.path.isfile(outfile) is True:
            file_suffix += 1
            outfile = file + '_' + str(file_suffix) + ext

        print "Warning: output file already exists." \
                " Outputting to new file", outfile

    if os.path.isfile(outfile) is True and newfile is False:
        pass
    else:
        ## Write proband labels as first line in output file. Only necessary
        ## when writing a new file.
        writefile(header_list, outfile)

    return outfile


## Currently not used
def appendcolumntofile(output, outfile):
    towrite = []
    for i in range(len(output)):
        towrite.append(map(str, output[i]))
    print towrite
    with open(outfile, 'r+') as outfile:
        for i, line in enumerate(outfile):
            columns = line.split("\t")
            print columns
            columns.insert(-1, towrite[i][0])
            columns.insert(-1, towrite[i][1])
            outfile.write("\t".join(columns) + "\n")


## Currently not used
def getlineage(ind, ind_dict, motherlist, fatherlist):
    lineage = [ind]
    indnum = ind_dict[ind]
    ind_father = fatherlist[indnum]
    ind_mother = motherlist[indnum]

    if ind_father != 0:
        lineage.extend(getlineage(ind_father, ind_dict, motherlist, fatherlist))

    if ind_mother != 0:
        lineage.extend(getlineage(ind_mother, ind_dict, motherlist, fatherlist))

    return lineage


## Currently not used
def getdescendants(ind, indlist, motherlist, fatherlist):
    descendants = [ind]

    if ind in motherlist and ind != 0:
        decnum = [i for i,n in enumerate(motherlist) if n == ind]
        for i in decnum:
            descendants.append(getdescendants(indlist[i], indlist, motherlist,
                                              fatherlist))
    elif ind in fatherlist and ind != 0:
        decnum = [i for i,n in enumerate(fatherlist) if n == ind]
        for i in decnum:
            descendants.append(getdescendants(indlist[i], indlist, motherlist,
                                              fatherlist))

    return descendants


def getprobands(indlist, motherdict, fatherdict):
    ## Probands are individuals who are not mothers or fathers
    probands = [indlist[i] for i,n in enumerate(indlist) if
                (indlist[i] not in motherdict and
                indlist[i] not in fatherdict)]

    indices = range(len(probands))
    proband_dict = dict(zip(probands, indices))

    return probands, proband_dict


def getleaves(indlist, motherlist, fatherlist):
    ## Leaves are those individuals without parents, indicated by a '0' in
    ## the pedigree.
    nodes = [indlist[i] for i,n in enumerate(motherlist) if n == 0]
    return nodes


def buildpath(indlist, ind_dict, offspring_dict, fatherlist, motherlist, startlist=None):
    ## Unless provided with a starting list, we start with the founders
    ## (leaves) of the pedigree
    if startlist is None:
        startlist = getleaves(indlist, motherlist, fatherlist)

    hasalleles = np.zeros(len(indlist)) # 0 - no alleles, 1 - has alleles
    path = [] # To store the individuals to receive alleles at each step
    nextgen = [] # Stores the individuals who can receive alleles next
    nextgen.extend(startlist) # First step is the founders

    ## Individuals for next step are chosen from offspring of previous step.
    ## Offspring with only one allele-calculated parent get carried over here.
    offspring_notselected = []


    while len(nextgen) > 0:
        ## Add step to path
        path.append(nextgen)

        ## Assign alleles to current step
        for ind in nextgen:
            hasalleles[ind_dict[ind]] = 1

        ## Find offspring of current step, to limit search size for next step,
        ## including offspring not calculated from previous step.
        offspring = []
        offspring.extend(offspring_notselected)
        offspring_notselected = []
        for ind in nextgen:
            try:
                offspring.extend(offspring_dict[ind])
            except KeyError: # When individual has no offspring
                continue
        offspring = list(OrderedDict.fromkeys(offspring))

        ## Find next step from offspring of current step.
        nextgen = []
        for off in offspring:
            ##@@ I think some of these checks are unnecessary
            if (hasalleles[ind_dict[off]] == 0 and
                motherlist[ind_dict[off]] != 0 and
                fatherlist[ind_dict[off]] != 0 and
                hasalleles[ind_dict[motherlist[ind_dict[off]]]] == 1 and
                    hasalleles[ind_dict[fatherlist[ind_dict[off]]]] == 1):
                nextgen.append(off)
            else:
                offspring_notselected.append(off)

    return path


def passalleles(indlist, ind_dict, fatherlist, motherlist, path, startlist=None):
    numremaining = len(indlist)
    ## Initialize alleles as empty lists. Start with leaves/ancestors, whose
    ## alleles we are tracking.
    alleles = [[] for i in range(len(indlist))]
    if startlist is None:
        startlist = getleaves(indlist, motherlist, fatherlist)

    numremaining -= len(startlist)

    ## Alleles for the leaves/ancestors are simply their own integer label.
    for ind in startlist:
        alleles[ind_dict[ind]] = [str(ind) + 'a', str(ind) + 'b']

    ## Each time we complete a loop, the alleles of more individuals can be
    ## computed. We continue until all individuals have alleles and the nextgen
    ## list is empty.
    for i in range(1, len(path)):
        numremaining -= len(path[i])
        for ind in path[i]:
            mother = motherlist[ind_dict[ind]]
            father = fatherlist[ind_dict[ind]]

            ## Each individual gets one allele from each parent, each passed
            ## down with a 50% chance.
            passallele1 = alleles[ind_dict[mother]][randint(2)]
            passallele2 = alleles[ind_dict[father]][randint(2)]
            alleles[ind_dict[ind]].extend([passallele1, passallele2])

    return zip(*alleles)


def run(datafile, outpath, Iterations, MakeNewFile, impsample_weights = None):
#    if MakeNewFile == 'overwrite':
#        outfile = outpath
#    else:
#        datafile_name, datafile_ext = os.path.splitext(os.path.basename(datafile))
#        outfile = os.path.join(outpath, datafile_name + "_out.csv")

    ## Prepares data for analysis
    datalist = loadpedfile(datafile)
    inds = datalist[0]
    fathers = datalist[1]
    mothers = datalist[2]
    ind_dict = datalist[3]
    father_dict = datalist[4]
    mother_dict = datalist[5]
    offspring_dict = datalist[6]

    ## Finds probands so we save their alleles
    probands, proband_dict = getprobands(inds, mother_dict, father_dict)
    proband_indices = [ind_dict[ind] for ind in probands]
    print len(probands), "probands out of", len(inds), "individuals in pedigree"

    ## Make list of alleles being passed
    allele_nums = getleaves(inds, mothers, fathers)
    alleles = []
    for allele in allele_nums:
        alleles.append(str(allele) + 'a')
        alleles.append(str(allele) + 'b')

    ## Check if output file already exists and create a new file if necessary.
    ##@@ Messy - if overwrite we need the path and filename, otherwise path only
    if MakeNewFile == 'overwrite':
        outfile = outpath
        newoutfile = outfile
        write_mode = 'w'
        writefile(probands, newoutfile, write_mode)
        write_mode = 'a'
    else:
        datafile_name, datafile_ext = os.path.splitext(os.path.basename(datafile))
        outfile = os.path.join(outpath, datafile_name + "_out.csv")
        newoutfile = createfile(outfile, probands, MakeNewFile)
        write_mode = 'a'

    ## Finds path to take through the pedigree
    print "Finding path through pedigree..."
    path = buildpath(inds, ind_dict, offspring_dict, fathers, mothers)
    print "Done"

    ## Runs simulations
    print "Running", Iterations, "iterations..."
    for i in range(Iterations):
        print "Iteration", i + 1, "/", Iterations
        ## Passes alleles from founders to probands
        results = passalleles(inds, ind_dict, fathers, mothers, path)

        ## We only want to save alleles of probands
        proband_allele0 = [results[0][i] for i in proband_indices]
        proband_allele1 = [results[1][i] for i in proband_indices]

        ## Writes output
        writefile(proband_allele0, newoutfile, write_mode)
        writefile(proband_allele1, newoutfile, write_mode)

    return newoutfile, alleles, probands, proband_dict, ind_dict, proband_allele0, proband_allele1
