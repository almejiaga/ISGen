import sys, os
import numpy as np
import pandas as pd
import tables
import argparse


def main(args):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--num-bins", metavar='|',
                        help="Number of evenly-spaced bins to use when " +\
                            "partitioning posterior probabilities",
                        type=int, default=11)

    ## DONE: Some of these are required only if file does not exist +p2 id:119
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-r", "--raw_liks_file", metavar='|',
                        help="File containing ancestor and likelihood for each" +\
                        " iteration within a simulation",
                        required=True)
    requiredNamed.add_argument("-o", "--outfile", metavar='|',
                        help="File to store output validation",
                        required=True)

    args = parser.parse_args()

    ## This is set here so that we can check whether we are calling from
    ## the CLI or another script above
    args.from_file = True
    args.CIs = None

    main(args)
