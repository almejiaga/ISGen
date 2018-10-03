import sys, os
sys.path.append('../bootstrap')
import validation_accuracy as validation
import itertools
import argparse
import tables
import pandas as pd
import bootstrap


def anc_from_panelfile(panel_file):
    """
    Returns the ancestor who generated the panel, by pulling
    out numbers from the file name. Robust when using validation
    pipeline.
    """
    path, filename = os.path.split(panel_file)
    anc = int(filter(str.isdigit, filename))

    return anc


def get_sym_anc(panel_file, symdat):
    """
    Returns the ancestor who generated the panel, after
    collapsing symmetries
    """
    gen_anc = anc_from_panelfile(panel_file)
    try:
        sym_gen_anc = int(symdat['Symmetrical'][symdat['Original']==gen_anc])
    except TypeError:
        sym_gen_anc = gen_anc

    return sym_gen_anc


def load_symmetries(symfile):
    """
    Loads symmetries into DataFrame
    """
    symdat = pd.DataFrame.from_csv(symfile)

    ## Split into original ancs and their symmetrical replacements
    symdat.columns = ['Symmetrical']
    symdat['Original'] = symdat.index
    symdat.reset_index(level=0, inplace=True)

    return symdat


def parse_paths(*path_files):
    """
    Returns list of file paths contained in CLI args
    """
    all_paths = []

    for path_file in path_files:
        paths = [line.strip() for line in open(path_file, 'r')]
        all_paths.append(map(os.path.expanduser, paths))

    return all_paths


def file_list_to_df(climb_files, control_files, sym_files):
    """
    Loads the likelihoods of all runs as a dict of DataFrames,
    and collapses symmetrical inds.
    """
    result_paths = zip(climb_files, control_files, sym_files)

    for i, (climb_file, control_file, sym_file) in enumerate(result_paths):
        try:
            B = bootstrap.Bootstrap(climb_file, control_file, sym_file)
        except tables.NoSuchNodeError:
            continue

        yield i, B.likdata


def write_hdf_node(key, df, outfile, sym_gen_anc=None):
    """
    Writes 'df' to node 'key' in hdf5 'outfile', with optional
    attribute labelling sym_gen_anc.
    """
    ## Connect to outfile and write array
    with pd.HDFStore(outfile) as store:
        label = 'run_' + str(key)
        store[label] = df
        if sym_gen_anc is not None:
            store.get_storer(label).attrs.sym_gen_anc = sym_gen_anc
    return


def main(args):
    ## Load data from CLI args
    result_paths = parse_paths(args.climb_files, args.control_files,
                               args.symmetry_files, args.panel_files)
    climb_files, control_files, sym_files, panel_files = result_paths
    data = file_list_to_df(climb_files, control_files, sym_files)

    ## Calculate the likelihood of each ancestor
    print "Loading likelihoods..."
    for i, df in data:
        print i
        ## Append indicator column denoting if ancestor generated panel
        symdat = load_symmetries(sym_files[i])
        sym_gen_anc = get_sym_anc(panel_files[i], symdat)

        ## Write results to file
        write_hdf_node(i, df, args.outfile, sym_gen_anc=sym_gen_anc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    ## Arguments for bootstrapping batch of likelihood files
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-o", "--outfile", metavar='|',
                        help="File to output total likelihoods for all runs, " +\
                        "as read from panel, climb, control, and symmetry " +\
                        "files", required=True)
    requiredNamed.add_argument("-p", "--panel-files", metavar='|',
                        help="File listing paths to panel files used for " +\
                        "validation", required=True)
    requiredNamed.add_argument("-l", "--climb-files", metavar='|',
                        help="File listing output files of allele " +\
                        "climbing simulations",
                        required=True)
    requiredNamed.add_argument("-c", "--control-files", metavar='|',
                        help="File listing control likelihoods for " +\
                         "the trees simulated in --climb-files",
                        required=True)
    requiredNamed.add_argument("-s", "--symmetry-files", metavar='|',
                    help="File listing symmetry files to use with " +\
                    "bootstrapping", required=True)

    args = parser.parse_args()

    main(args)
