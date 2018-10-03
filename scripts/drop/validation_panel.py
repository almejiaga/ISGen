# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:37:36 2015

@author: dominic
"""
import numpy as np
import sys, os
from collections import defaultdict, Counter
import argparse
import pandas as pd
import attr


def build_region_dict(regionfile, all_probands, delimiter='\t'):
    """ Returns a dict with a set of probands belonging to each region """
    header = ['Proband', 'Region']
    regions = pd.read_csv(regionfile, sep=delimiter, names=header,
                          index_col=0, dtype=dict(zip(header, [int, str])))

    region_dict = {}
    for region, probands in regions.groupby('Region').groups.iteritems():
        region_dict[region.strip()] = set(probands).intersection(all_probands)

    # Add a region containing all probands
    region_dict['All Probands'] = set(regions.index).intersection(all_probands)

    # Assign a region of 'None' if the proband is not listed in the file
    no_region = set(all_probands).difference(region_dict['All Probands'])
    if len(no_region) > 0:
        region_dict['None'] = no_region

    return region_dict


def load_simfile(probands, simfile, startsim, num_sims):
    """
    Generator yielding a single simulation results on each call, each of which
    contains two alleles per individual.
    """
    # Sims are two lines each, plus header
    startline = np.max([(startsim * 2) - 1, 1])
    df = pd.read_csv(simfile, index_col=False, chunksize=2, header=None,
                     skiprows=startline, names=probands, dtype=str)

    # Return the transpose of the simulation data, so that we can group
    # inds who share alleles more easily
    i = 0
    for sim in df:
        i += 1
        yield sim.transpose()

        if i == num_sims:
            break


def group_alleles(sim_data):
    """
    Returns groups of individuals who share alleles
    """
    allele_groups = defaultdict(Counter)

    # Group inds who share alleles
    allele_labels = list(sim_data)
    assert len(allele_labels) == 2

    allele0 = sim_data.groupby(allele_labels[0])
    allele1 = sim_data.groupby(allele_labels[1])

    # Combine alleles of each ind
    for key, value in allele0.groups.iteritems():
        allele_groups[key] += Counter(value)

    for key, value in allele1.groups.iteritems():
        allele_groups[key] += Counter(value)

    return allele_groups


def check_attribs(instance, attribute, value):
    try:
        if attribute.name == 'hom_het':
            assert value in ['hom', 'het']
    except:
        print "Validation error!"
        import pdb; pdb.set_trace()
        raise


@attr.s
class sorted_panels(object):
    """
    Returns panels within the specified size range, along with the regional
    allele frequencies (given by the individuals not included in the panel).
    """
    # TODO: Allow different sampling strategies, ie per-region +p3 id:149
    allele_groups = attr.ib()
    region_dict = attr.ib()
    regions = attr.ib()
    min_size = attr.ib()
    max_size = attr.ib()

    sample_size = attr.ib(default=1000)
    first_panel = attr.ib(default=False)
    hom_het = attr.ib(default='hom', validator=check_attribs)
    min_carriers = attr.ib(default=1)

    ancs = attr.ib(init=False, default=attr.Factory(list))
    panels = attr.ib(init=False, default=attr.Factory(list))
    regional_counts = attr.ib(init=False, default=attr.Factory(list))

    def __attrs_post_init__(self):
        all_inds = list(self.region_dict['All Probands'])

        for allele, ind_counts in self.allele_groups.iteritems():
            # Pull inds out of Counter object
            inds = ind_counts.keys()

            # Only create one panel per anc if flag is set
            # TODO: Does nothing when processing one sim at a time +p3 id:171
            anc = allele[:-1]
            if anc in set(self.ancs):
                continue

            # NOTE: We assume we find all homozygotes in the population id:173
            # TODO: Implement sampling strategy here as well +p3 id:172
            if self.hom_het == 'hom':
                panel = [ind for ind, count in ind_counts.iteritems()
                         if count == 2]
                not_panel = set(inds).difference(panel)
            else:
                # Randomly choose a panel of known carriers from all probands
                # NOTE: This considers homs as hets +p3 id:175
                self.sample_size = np.min([int(np.round(len(inds) / 2.)),
                                            np.random.randint(self.min_size,
                                                        self.max_size + 1)])
                panel = np.random.choice(inds, size=self.sample_size,
                                          replace=False)
                not_panel = set(inds).difference(panel)

            # Only include groups of proper size
            if len(panel) < self.min_size or len(panel) > self.max_size:
                continue

            # Make sure we have the minimum number of non-panel carriers
            if len(not_panel) < self.min_carriers:
                continue

            # Calculate regional allele frequencies, including total number
            # of inds in each region
            counts = []
            for r in self.regions:
                count = len(self.region_dict[r].intersection(not_panel))
                ninds_region = len(self.region_dict[r])
                counts_to_write = str(count) + '/' + str(ninds_region)
                counts.append(counts_to_write)

            # Append results
            self.ancs.append(anc)
            self.panels.append(panel)
            self.regional_counts.append(counts)


def write_panels(ancs, panels, regional_counts, regions, outfile):
    """
    Write results to file in format:

        gen_anc \t clause \t region1_freq \t ... \t regionN_freq \n
    """
    # NOTE: This would be more flexible as HDF5 file +p3 id:174
    with open(outfile, 'w') as f:
        # Write header
        header = 'Generating Anc' + '\t' + 'Clause' + '\t'
        for region in regions:
            header += str(region) + '\t'

        # Replace last tab with newline and write to file
        assert header[-1] == '\t'
        header = header[:-1] + '\n'
        f.write(header)

        for anc, panel, counts in zip(ancs, panels, regional_counts):
            # Add gen_anc
            towrite = str(anc) + '\t'

            # Add panel, which is comma-delimited
            towrite += ','.join(map(str, panel)) + '\t'

            # Add region allele counts
            assert len(counts) == len(regions)
            for count in counts:
                towrite += str(count) + '\t'

            # Change last tab of allele counts to a newline
            assert towrite[-1] == '\t'
            towrite = towrite[:-1] + '\n'

            # Write to file
            f.write(towrite)


def main(args):
    # Parse CLI arguments
    simfile = os.path.expanduser(args.simfile)
    outfile = os.path.expanduser(args.outfile)

    if args.delimiter is None:
        delimiter = '\t'
    else:
        delimiter = args.delimiter

    if args.regionfile is not None:
        regionfile = os.path.expanduser(args.regionfile)
    else:
        regionfile = None

    if args.panel_ancs_file is not None:
        ancs_for_panel = map(int, np.genfromtxt(args.panel_ancs_file,
                                                delimiter=',',
                                                skip_header=False))
    else:
        ancs_for_panel = []

    panel_minsize, panel_maxsize = map(int, args.minmax.split(','))

    # Probands are stored as the header of the file
    probands = np.genfromtxt(simfile, max_rows=1, dtype=int, delimiter=',')

    # Load regions from regionfile
    print "Loading regions..."
    if regionfile is not None:
        region_dict = build_region_dict(regionfile, probands, delimiter)
    else:
        region_dict = {'All Probands': set(probands)}
    regions = sorted(region_dict.keys())
    print len(region_dict), "regions found"

    # Load simulation results as a generator
    sims = load_simfile(probands, simfile, args.startsim, args.num_sims)

    # Iterate through simulations in file
    ancs, panels, regional_counts = [], [], []
    for sim_data in sims:
        # Sort simulation results into clauses
        print "Finding panels..."
        allele_groups = group_alleles(sim_data)

        # Separate groups of allele-sharing inds into panel (connected to
        # genealogy) and non-panel (disconnected) inds
        parsed_panels = sorted_panels(allele_groups, region_dict, regions,
                                      panel_minsize, panel_maxsize,
                                      hom_het=args.hom_het,
                                      min_carriers=args.min_carriers)

        # Append results
        ancs.extend(parsed_panels.ancs)
        panels.extend(parsed_panels.panels)
        regional_counts.extend(parsed_panels.regional_counts)

    if len(panels) == 0:
        print "No clauses with length in range", panel_minsize, "to", +\
                                                        panel_maxsize
        sys.exit()

    # Write panels to outfile. First column is the generating ancestor,
    # followed by a comma-delimited list of the probands who received their
    # alleles, followed by regional allele frequencies
    print "Writing panels to file..."
    write_panels(ancs, panels, regional_counts, regions, outfile)
    print "Done!"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--panel-ancs-file", metavar='|',
                        help="File containing ancestors to create panels for")
    parser.add_argument("-i", "--first-panel",
                        help="Only store the first panel found for each " +\
                        "ancestor",
                        action='store_true')
    parser.add_argument("-r", "--regionfile", metavar='|',
                        help="File listing the region of each proband, " +\
                       "with no header")
    parser.add_argument("-d", "--delimiter", metavar='|',
                        help="Delimiter used in --regionfile. If space, " +\
                        "surround by double quotes. Default is tab.")
    parser.add_argument("-c", "--min-carriers", metavar='|',
                        help="Minimum number of non-panel inds required" +\
                        " to calculate allele frequency",
                        type=int, default=1)
    parser.add_argument("-s", "--startsim", metavar='|',
                        help="Simulation number in --simfile from which " +\
                        "to start searching for panels",
                        type=int, default=0)
    parser.add_argument("-n", "--num-sims", metavar='|',
                        help="Number of simulations in --simfile to " +\
                        "process",
                        type=int, default=-1)
    # DONE: Reimplement creation of hom/het panels +p3 id:151
    parser.add_argument("-H", "--hom-het", metavar='|',
                        help="Indicates whether hom or het panels will be " +\
                        "constructed. One of ['hom' (default), 'het']",
                        choices=['hom', 'het'], default='hom')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-f", "--simfile", metavar='|',
                        help="File containing results of allele dropping " +\
                        "simulations.",
                        required=True)
    requiredNamed.add_argument("-o", "--outfile", metavar='|',
                        help="File to store output, including panels, " +\
                        "generating ancestors, and regional allele frequencies",
                        required=True)
    requiredNamed.add_argument("-m", "--minmax", metavar='|',
                        help="Minimum and Maximum size of panels to create" +\
                            ", in form min,max (no spaces)",
                        required=True)


    args = parser.parse_args()

    main(args)
