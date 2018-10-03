import sys, os
import traceback
import glob
import itertools
import subprocess32
import numpy as np
import tables
import pandas as pd
import attr
import argparse


@attr.s
class H5Tools(object):
    """
    Loads the data and attributes of a simple hdf5 file. Note that
    the structure must consist only of leaves descended from the
    root group for copying to be complete.
    """
    h5file = attr.ib()
    global_attrs = attr.ib(init=False)
    load_data = attr.ib(default=True)
    load_posteriors = attr.ib(default=False)
    arrays = attr.ib(init=False, default=attr.Factory(list))
    atoms = attr.ib(init=False, default=attr.Factory(list))
    array_attrs = attr.ib(init=False, default=attr.Factory(list))
    filters = attr.ib(init=False, default=attr.Factory(list))
    names = attr.ib(init=False, default=attr.Factory(list))
    liks_only = attr.ib(default=False)

    def __attrs_post_init__(self):
        with tables.open_file(self.h5file, 'r') as f:
            ## Store the attributes of the root group itself
            self.global_attrs = f.root._v_attrs

            ## Store the data and attributes of each Leaf within
            ## the group
            for node in f.walk_nodes(f.root, classname='Leaf'):
                ## If flag is set, remove posteriors since they can't
                ## be simply appended
                if not self.load_posteriors:
                    if node.name == 'posteriors':
                        continue

                ## Load only 'ancs' and 'liks' if flag is set
                if self.liks_only and node.name not in ['ancs', 'liks', 'tot_liks',
                                                        'control_liks']:
                    continue

                ## Store array attributes
                self.array_attrs.append(node._v_attrs)
                self.names.append(node.name)
                self.filters.append(node.filters)
                self.atoms.append(node.atom)

                ## Load data if flag is set
                if self.load_data:
                    self.arrays.append(node[:])
                else:
                    self.arrays.append([])


    def subsample(self, nrows):
        """
        Draws a random sample of rows from the file
        """
        total_rows = self.arrays[0].shape[0]
        rows = np.random.choice(np.arange(total_rows), size=nrows,
                                replace=False)

        sub_sample = []
        for array in self.arrays:
            try:
                sub_sample.append(array[rows])
            except TypeError:
                ## If we can't specify rows as a slice, use
                ## list comprehension
                sub_sample.append([array[i] for i in rows])

        return sub_sample


    def get_dtype(self, array):
        """
        Returns the dtype of the given array, which will either be
        a normal numpy array, a numpy array with 'object' datatype, or
        a list of numpy arrays
        """
        if type(array) == list:
            dtype = array[0].dtype
        elif array.dtype == 'object':
            dtype = array[0].dtype
        else:
            dtype = array.dtype

        return dtype


    def add_array(self, f, array_index, data):
        """
        Creates an array of the specified type in root node of the open file
        object provided
        """
        name = self.names[array_index]
        array_class = self.array_attrs[array_index].CLASS
        filters = self.filters[array_index]

        ## Create the proper atom for the array, based on dtype
        dtype = self.get_dtype(data)
        atom = tables.Atom.from_dtype(dtype)

        ## A vlarray is easiest to create line-by-line
        if array_class == "VLARRAY":
            vl = f.create_vlarray(f.root, name, atom=atom, filters=filters)
            for row in data:
                vl.append(row)
        else:
            f.create_array(f.root, name, obj=data, atom=atom)


    def write(self, out_file, array_data):
        with tables.open_file(out_file, 'w') as f:
            for i in range(len(self.arrays)):
                self.add_array(f, i, array_data[i])


@attr.s
class H5Merge(object):
    """
    Tools for merging a list of hdf5 files
    """
    h5_file_list = attr.ib()
    liks_only = attr.ib(default=False)


    def load_files(self):
        for h5_file in self.h5_file_list:
            file_contents = H5Tools(h5_file, liks_only=self.liks_only)

            yield file_contents


    def merge_files(self, loaded_file1, loaded_file2):
        """
        Merges the second provided file into the first one
        """
        for i in range(len(loaded_file1.arrays)):
            assert loaded_file1.names[i] == loaded_file2.names[i]
            
            ## Load and stack arrays
            array1 = loaded_file1.arrays[i]
            array2 = loaded_file2.arrays[i]
            loaded_file1.arrays[i] = np.append(array1, array2, axis=0)

        return loaded_file1


    def merge_all(self):
        ## Create an empty file by loading the first file without the data
        merged_file = H5Tools(self.h5_file_list[0], load_data=False,
                              liks_only=self.liks_only)

        ## Load all other files and merge them
        for loaded_file in self.load_files():
            self.merge_files(merged_file, loaded_file)

        return merged_file


def find_files_to_merge(path):
    """
    Finds simulation files to merge within the provided path
    """
    control_file = climb_file = None

    for f in glob.glob(path + "/*"):
        if "control" in f:
            control_file = f
        elif "climb" in f:
            climb_file = f

    assert control_file is not None and climb_file is not None

    return climb_file, control_file


def copy_sim_params(path, out_path):
    """
    Copies files used in the simulation present in 'path' to 'outpath'
    """
    out_path = os.path.expanduser(out_path)
    assert os.path.exists(out_path)
    assert os.path.exists(path)

    command = 'cp -Rv ' + path + '/*.txt ' + out_path + '/'
    subprocess32.call([command], shell=True)

    return


def make_if_not_exists(path):
    """
    Creates 'path' unless it exists, in which case do nothing.
    """
    try:
        os.makedirs(path)
    except OSError, e:
        ## If error raised because path exists, do nothing
        if e.errno != 17:
            raise


def match_merged_output(paths_to_merge, merged_paths):
    """
    Returns ordered lists of matching simulation files to merge, along with
    the full file path to output the merged file
    """
    ##HACK: List of empty lists not very elegant +p3 id:202
    files_to_merge = [ [] for i in range(len(paths_to_merge))]
    new_merged_paths = []

    for file_params in zip(merged_paths, *paths_to_merge):
        merged_path, to_merge = file_params[0], file_params[1:]
        new_paths = []

        ## Create merged path if it doesn't exist
        make_if_not_exists(merged_path)

        for sim_path, l in zip(to_merge, files_to_merge):
            ## Add to list of path lists, and add matching output paths
            climb_file, control_file = find_files_to_merge(sim_path)
            
            ## Append the results to lists in files_to_merge
            l.extend([climb_file, control_file])

            ## Copy other files in the directory to the new path
            copy_sim_params(sim_path, merged_path)

        ## Add appropriate output file for each set of files to be merged
        out_fname = lambda x: os.path.join(merged_path, os.path.split(x)[1])
        new_paths = [out_fname(f) for f in [climb_file, control_file]]
        new_merged_paths.extend(new_paths)

    ## We should have the same number of merged files and output dirs
    assert set(map(len, files_to_merge)) == set([len(new_merged_paths)])

    return files_to_merge, new_merged_paths


def main(args):
    ## Read list of paths to write merged files to
    merged_list = os.path.expanduser(args.outfiles)
    with open(merged_list, 'r') as f:
        l = [line.strip() for line in f]
        all_merged_paths = map(os.path.expanduser, l)

    ## Read list of paths in each file
    path_lists = map(os.path.expanduser, args.path_lists)
    paths_to_merge = []
    for path_list in path_lists:
        with open(path_list, 'r') as f:
            l = [line.strip() for line in f]
            paths_to_merge.append(map(os.path.expanduser, l))

    ## Find files to merge within provided simulation paths, and modify
    ## output file list so that multiple files can be written to the same
    ## directory
    files_to_merge, merged_files = match_merged_output(paths_to_merge,
                                                       all_merged_paths)

    ## Iterate through file to merge and their associated output files
    for file_params in zip(merged_files, *files_to_merge):
        try:
            merged_file, to_merge = file_params[0], file_params[1:]

            print "Merging", to_merge, "into", merged_file
            H = H5Merge(to_merge, liks_only=args.liks_only)
            m = H.merge_all()

            ## Subsample merged data if n_rows is specified
            if args.n_rows > 0:
                data = m.subsample(args.n_rows)
            else:
                data = m.arrays

            ## Write merged file
            print "Writing merged output"
            m.write(merged_file, data)

        ## Log errors and continue without killing script
        except KeyboardInterrupt:
            raise
        except:
            print "Error!"
            with open('err.log', 'a') as f:
                f.write('-'*60)
                f.write('Error processing the following files:\n')
                f.write(str(merged_file) + '\n')
                for _file in to_merge:
                    f.write(str(_file) + '\n')
                traceback.print_exc(file=f)
            continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--n_rows', metavar='',
                        help='Number of rows to sample (without replacement) from' +\
                        ' the merged file', type=int, default=-1)
    parser.add_argument('-p', '--path_lists', metavar='', nargs='+',
                        help='Paths to files which contain ordered lists of files ' +\
                        'to be merged')
    parser.add_argument('-o', '--outfiles', metavar='',
                        help='File containing ordered list of output paths for ' +\
                        'the merged files')
    parser.add_argument('-l', '--liks_only',
                        help='Flag to merge only ancestors and likelihoods, ' +\
                        'without corresponding trees and genotypes, with ' +\
                        'better performance', action='store_true')

    args = parser.parse_args()

    main(args)
