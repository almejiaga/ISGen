import numpy as np
import sys,os
sys.path.append('../climb/')
import ped
import argparse
import pandas as pd

class PedSymmetry:
    def __init__(self, P):
        self.P = P


    def get_parents(self, ind):
        mother = self.P.mothers[self.P.ind_dict[ind]]
        father = self.P.fathers[self.P.ind_dict[ind]]

        return mother, father


    def offspring_root_anc(self, ind, commonancs, all_lineages):
        have_path = []
        for off in self.P.offspring_dict[ind]:
            if off in all_lineages:
                have_path.append(off)

        offspring_root = []
        if len(have_path) == 1 and have_path[0] in commonancs:
            offspring_root = have_path

        return offspring_root


    def isfounder(self, ind):
        founder = 0
        mother, father = self.get_parents(ind)
        if mother == 0:
            founder += 0.5
        if father == 0:
            founder += 0.5

        return founder


    def get_commonancs(self, probands):
        all_lineages = set(self.P.ordered_lineage(probands[0]).keys())
        commonancs = set(self.P.ordered_lineage(probands[0]).keys())

        for proband in probands[1:]:
            lineage = set(self.P.ordered_lineage(proband).keys())
            all_lineages = all_lineages.union(lineage)
            commonancs = commonancs.intersection(lineage)

        return list(commonancs), all_lineages


    def find_partners(self, ind):
        children = self.P.offspring_dict[ind]
        partners = set()
        for child in children:
            partners.update(self.get_parents(child))

        return tuple(sorted(set(partners)))


    def find_monogamous_partners(self, ind):
        partners = self.find_partners(ind)
        partner_partners = set()
        for partner in partners:
            partner_partners.add



    def find_monogamous_in_group(self, group):
        monogamous = set()
        group = set(group)
        for ind in group:
            partners = self.find_partners(ind)
            if len(partners) == 2 and len(group.intersection(partners)) == 2:
                monogamous.add(ind)

        return list(monogamous)


    def get_root_ancs(self, probands):
        commonancs, all_lineages = self.get_commonancs(probands)
        current_root_ancs = [x for x in commonancs if self.isfounder(x) == 1]

        root_ancs = set()
        while True:
            new_root_ancs = []
            to_remove = []
            for f in current_root_ancs:
                root_offspring = self.offspring_root_anc(f, commonancs, all_lineages)

                ## If root has offpsring who is also a root, we take the offspring
                ## and remove the parent
                if len(root_offspring) > 0:
                    to_remove.append(f)

                new_root_ancs.extend(root_offspring)

            if len(new_root_ancs) == 0:
                break

            root_ancs.update(new_root_ancs)
            current_root_ancs = new_root_ancs
            root_ancs = root_ancs.difference(to_remove)

        return root_ancs


    def get_symm_inds(self, probands):
        root_ancs = get_root_ancs(probands)

        sym_dict = {}
        for anc in root_ancs:
            for sym in self.P.ordered_lineage(anc).keys():
                sym_dict[sym] = anc

        return sym_dict


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
    P = ped.Pedigree(pedfile)
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
