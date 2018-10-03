import numpy as np
import sys,os
sys.path.append(os.path.expanduser('~/project/GENLIB/bin/scripts/anc_finder'))
import anc_finder_ped as ped
import argparse
import pandas as pd

class PedSymmetry:
    def __init__(self, P):
        self.P = P


    def get_parents(self, ind):
        mother = self.P.mothers[self.P.ind_dict[ind]]
        father = self.P.fathers[self.P.ind_dict[ind]]

        return mother, father


    def ascending_sibs(self, ind, commonancs):
        mother, father = self.get_parents(ind)
        Msibs = set(self.P.offspring_dict[mother])
        Fsibs = set(self.P.offspring_dict[father])

        if mother == 0:
            Msibs = set([ind])
        if father == 0:
            Fsibs = set([ind])

        sibs = Msibs.union(Fsibs)

        return sibs.intersection(set(commonancs))


    def isfounder(self, ind):
        founder = 0
        mother, father = self.get_parents(ind)
        if mother == 0:
            founder += 0.5
        if father == 0:
            founder += 0.5

        return founder


    def get_commonancs(self, probands):
        commonancs = set(self.P.ordered_lineage(probands[0]).keys())
        for proband in probands[1:]:
            lineage = set(self.P.ordered_lineage(proband).keys())
            commonancs = commonancs.intersection(lineage)

        return list(commonancs)


    def find_partners(self, ind):
        children = self.P.offspring_dict[ind]
        partners = []
        for child in children:
            partners.extend(self.get_parents(child))

        return tuple(sorted(set(partners)))


    def find_monogamous(self, founders):
        is_monogamous = lambda x: len(x) == 2
        monogamous = filter(is_monogamous, map(self.find_partners, founders))

        return list(set(monogamous))


    def get_symm_inds(self, probands, combine=True):
        commonancs = self.get_commonancs(probands)

        not_founders = set(filter(lambda x: self.isfounder(x) == 0, commonancs))
        only_child = lambda x: len(self.ascending_sibs(x, commonancs)) == 1
        no_sibs = set(filter(only_child, commonancs))

        root_ancs = not_founders.intersection(no_sibs)

        sym_dict = {}
        for anc in root_ancs:
            for sym in self.P.ordered_lineage(anc).keys():
                sym_dict[sym] = anc

        founders = set(filter(lambda x: self.isfounder(x) != 0, commonancs))
        ## If they're monogamous and one of them is symmetrical, the other must
        ## be as well
        more_than_one_child = lambda (x,y): x not in sym_dict
        monogamous = filter(more_than_one_child, self.find_monogamous(founders))

        print len(commonancs), "common ancestors"
        print len(sym_dict), "inds with some symmetry"
        print len(monogamous), "monogamous founder couples"

        ## Add monogomous couples to symmetry dict, and return it
        ## if specified
        if combine is True:
            for x, y in monogamous:
                ## Need to match which monogamous ind is the target, if
                ## they are already in the symmetry dict. 
                if x in sym_dict:
                    sym_dict[x] = y
                else:
                    sym_dict[y] = x

            return sym_dict
        else:
            return sym_dict, monogamous


    def write_symmetrical_inds(self, probands, outfile):
        sym_dict = self.get_symm_inds(probands, combine=True)
        if len(sym_dict) > 0:
            sym_inds = pd.DataFrame.from_dict(sym_dict, orient='index')
            sym_inds.columns = [self.P.pedfile]
            sym_inds.index.name = probands
            sym_inds.to_csv(outfile)

            return 0

        return 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile",
                        help="Output file. Both symmetrical inds and " +\
                        "monogamous couples will be saved in the same file")
    parser.add_argument("-s", "--split",
                        help="Display monogamous couples separately from " +\
                        "symmetrical individuals. Any file written will " +\
                        "still combine the two")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-f", "--pedfile",
                        help="File containing the pedigree to analyze",
                        required=True)

    requireOne = parser.add_argument_group('require one of')
    group = requireOne.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--probandfile",
                        help="File containing column of probands")
    group.add_argument("-l", "--probandlist",
                        help="List of probands separated by commas (no spaces)")
    group.add_argument("-P", "--probandpaths-file",
                        help="File containing paths of multiple proband files." +\
                                " Writes symmetries of each panel to the same " +\
                                "directory, in a file named 'symmetry.txt'")

    args = parser.parse_args()

    if args.probandpaths_file and args.outfile:
        print "Incompatible command-like options: Cannot specify outfile for", +\
                "batch symmetry output."
        sys.exit()

    pedfile = os.path.expanduser(args.pedfile)

    if args.probandfile is not None:
        probands = map(int,[line.strip() for line in open(args.probandfile, 'rU')])
    if args.probandlist is not None:
        probands = map(int, args.probandlist.split(','))
    if args.probandpaths_file is not None:
        allprobands = []
        outpaths = []
        probandpaths = [line.strip() for line in
                                        open(args.probandpaths_file, 'rU')]
        for path in probandpaths:
            allprobands.append(map(int, [line.strip() for line in
                                            open(path, 'rU')]))
            outpaths.append(os.path.split(path)[0] + '/symmetry.txt')


    print "Loading pedigree from file:", pedfile
    P = ped.anc_finder_pedigree(pedfile)
    print "Done."
    print "Finding symmetries..."
    PS = PedSymmetry(P)

    if args.outfile is not None:
        print "Probands:", probands
        print "Outputting symmetries to", args.outfile
        PS.write_symmetrical_inds(probands, args.outfile)
    elif args.probandpaths_file is not None:
        for probands, path in zip(allprobands, outpaths):
            print "Finding symmetries for probands:", probands
            PS.write_symmetrical_inds(probands, path)
            print "Symmetries written to:", path
    else:
        print "Probands:", probands
        print PS.get_symm_inds(probands)


if __name__ == "__main__":
    main()
