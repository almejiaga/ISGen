import sys, os
import subprocess32
import numpy as np
import argparse


def batch_gen_anc(batch_dir):
    """
    Reads the generating ancestor from the name of the provided panel file
    """
    panel_file = [f for f in os.listdir(batch_dir) if 'panel' in f][0]
    anc = int(panel_file.split('_')[0])

    return anc


def get_batch_freq_file(batch_dir):
    """
    Returns the file containing regional allele frequencies associated with
    the simulated panel
    """
    freq_file = [f for f in os.listdir(batch_dir) if 'freq' in f][0]

    return freq_file


def load_expected_freqs(expected_freqs_file):
    """
    Returns a list of the expected allele frequencies for each region
    """
    freqs = []
    with open(expected_freqs_file, 'r') as f:
        for line in f:
            ## Recover freq from string formatted as
            ## 'Region: num_alleles/ninds_region'
            region_freq = line.split(' ')[-1].split('/')[0]
            freqs.append(region_freq)

    return freqs


def same_regional_freqs(regional_freqs_list):
    """
    Compares the entries of each list of regional allele freqs contained
    in regional_freqs_list. Returns True only if all are identical.
    """
    for region_freqs in zip(*regional_freqs_list):
        if len(set(region_freqs)) != 1:
            return False

    return True


def main(args):
    ## Read output paths to be compared
    out_paths1 = [line.strip() for line in open(args.out_paths_file1, 'r')]
    out_paths2 = [line.strip() for line in open(args.out_paths_file2, 'r')]
    out_paths1 = map(os.path.expanduser, out_paths1) 
    out_paths2 = map(os.path.expanduser, out_paths2) 

    ## To store jobs which have the same generating ancestor and regional
    ## allele frequencies
    matching_jobs = []

    for i, p1 in enumerate(out_paths1):
        ## Filter out higher-level dirs
        ##TODO: This and same below break with trailing forward slash +p2 id:205
        if "_" not in p1.split('/')[-1]:
            continue
        
        with open(os.devnull, 'w') as FNULL:
            for j, p2 in enumerate(out_paths2):
                ## Filter out higher-level dirs
                if "_" not in p2.split('/')[-1]:
                    continue

                ## Find the ancestor who generated the panel for each batch job
                anc1 = batch_gen_anc(p1)
                anc2 = batch_gen_anc(p2)

                freqs_file1 = os.path.join(p1, get_batch_freq_file(p1))
                freqs_file2 = os.path.join(p2, get_batch_freq_file(p2))
                freqs1 = load_expected_freqs(freqs_file1)
                freqs2 = load_expected_freqs(freqs_file2)

                if not same_regional_freqs([freqs1, freqs2]):
                    continue
                
                matching_jobs.append((p1, p2))
                print len(matching_jobs), "matches out of", i, "jobs"

                with open(os.path.join(args.outpath, 'jobs1.txt'), 'a') as f:
                    f.write(p1 + '\n')
                with open(os.path.join(args.outpath, 'jobs2.txt'), 'a') as f:
                    f.write(p2 + '\n')

        ## Here we know that p2=out_paths2[j] has been added, so we can remove it
        if args.unique:
            del out_paths2[j]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('out_paths_file1',
                        help='File listing output paths of first batch' +\
                        ' of simulations')
    parser.add_argument('out_paths_file2',
                        help='File listing output paths of second batch' +\
                        ' of simulations')
    parser.add_argument('outpath',
                        help='Path in which to write files containing paths ' +\
                        'of matching jobs, in two sorted text files, "jobs1.txt" ' +\
                        'and "jobs2.txt"')
    parser.add_argument('-u', '--unique',
                        help='Unique merge: Each path will be merged only once',
                        action='store_true')

    args = parser.parse_args()

    main(args)
