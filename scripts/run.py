import sys, os
sys.path.append('control/')
sys.path.append('climb/')
sys.path.append('msub/')
import argparse, ConfigParser
import coalesce, ind_contrib, control_lik, regional_freq

## Older versions of scipy.stats.hypergeom.pmf return nan instead of 0
## in certain cases - affects control likelihoods
## HACK: This might not be a very robust way of version checking id:176
import scipy as sp
print "Scipy version:", sp.__version__
major, minor, revision = [int(x) for x in sp.__version__.split('.')]
if major == 0 and minor < 17:
    assert sp_version >= 17, "Minimum scipy 0.17, found " + sp.__version__

def main(args):
    ## Run contribution allele dropping simulations
    if not os.path.isfile(os.path.expanduser(args.contrib_outfile)):
        assert args.contrib_iterations is not None, "Missing iterations " +\
                "for allele dropping simulations"

        contrib_args = argparse.Namespace(
                regionfile=args.proband_regions,
                verbose=args.verbose,
                pedfile=args.pedfile,
                outfile=args.contrib_outfile,
                iterations=args.contrib_iterations)

        ind_contrib.main(contrib_args)

    ## Run coalescence simulations if file does not already exist
    if not os.path.isfile(os.path.expanduser(args.climb_outfile)):
        assert args.climb_iterations is not None, "Missing iterations " +\
                "for allele climbing simulations"

        coalesce_args = argparse.Namespace(
                sim_homs=args.sim_homs,
                pedfile=args.pedfile,
                outfile=args.climb_outfile,
                iterations=args.climb_iterations,
                kinship=args.kinship,
                verbose=args.verbose)

        if 'hetfile' in vars(args):
            vars(coalesce_args)['hetfile'] = args.hetfile
        if 'homfile' in vars(args):
            vars(coalesce_args)['homfile'] = args.homfile

        coalesce.main(coalesce_args)

    ## Calculate control likelihoods of each simulated tree
    if (not os.path.isfile(os.path.expanduser(args.control_outfile)) and
            os.path.isfile(os.path.expanduser(args.contrib_outfile))):

        control_args = argparse.Namespace(
                pedfile=args.pedfile,
                climb_likfile=args.climb_outfile,
                outfile=args.control_outfile,
                contrib_file=args.contrib_outfile,
                freq_params=args.freq_params,
                verbose=args.verbose)

        control_lik.main(control_args)

    ## Calculate estimated regional allele frequencies
    if (args.regional_freq_outfile is not None
        and (not os.path.isfile(os.path.expanduser(args.regional_freq_outfile))
             and os.path.isfile(os.path.expanduser(args.contrib_outfile)))):

        N, n, k = [float(x) for x in args.freq_params.split(',')]
        observed_freq = 1. * k / n

        regional_freq_args = argparse.Namespace(
                regionsfile=args.regions_to_calculate,
                pedfile=args.pedfile,
                contrib_h5=args.contrib_outfile,
                climb_likfile=args.climb_outfile,
                observed_freq=observed_freq,
                control_freq_likfile=args.control_outfile,
                freq_params=args.freq_params,
                outfile=args.regional_freq_outfile)

        regional_freq.main(regional_freq_args)


if __name__ == "__main__":
    ## Climbing likelihood args
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--label", metavar='|',
                        help="Label to be prepended to output files")
    parser.add_argument("-s", "--sim_homs", metavar='|', default=0.5,
                        type=float, help="Frequency with which hets become " +\
                        "homs when climbed to. (default: 0.5)")
    parser.add_argument("-P", "--proband_regions", metavar='|',
                        help="File listing the region of each proband, " +\
                            "with no header. By default all probands are " +\
                            "assumed to be within a single region")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-r", "--regions_to_calculate", metavar='|',
                        help="""File listing regions for which to calculate
                                expected allele frequencies""")
    parser.add_argument("-i", "--contrib_iterations", metavar='|',
                        help="Number of allele drops to perform from each " +\
                            "individual in the pedigree",
                        type=int)
    parser.add_argument("-I", "--climb_iterations", metavar='|',
                        help="Number of allele climbing iterations to " +\
                            "perform from the provided homs and hets",
                        type=int)
    parser.add_argument("-k", "--kinship",
                        action='store_true',
                        help="Flag for use of kinship-based importance sampling,\n" +\
                             " which trades per-iteration computation time for\n" +\
                             " faster per-iteration convergence. Performance \n" +\
                             "increase only seen downstream, for control \n" +\
                             "likelihoods or regional allele frequency estimates")

    ## DONE: Some of these are required only if file does not exist +p2 id:119
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-b", "--contrib_outfile", metavar='|',
                        help="File to store allele dropping simulations " +\
                        " output in hdf5 format",
                        required=True)
    requiredNamed.add_argument("-p", "--pedfile", metavar='|',
                        help="Pedigree file",
                        required=True)

    single_required = parser.add_argument_group('Required if running on ' +\
                                                'single panel')
    single_required.add_argument("-o", "--climb_outfile", metavar='|',
                        help="Directory in which to store climbing output")
    single_required.add_argument("-n", "--control_outfile", metavar='|',
                        help="File to store control likelihoods for each " +\
                        "tree, in hdf5 format")
    single_required.add_argument("-f", "--freq_params", metavar='|',
                        help="Parameters describing observations of the " +\
                            "population allele frequency, in the format " +\
                            "'pop_size,num_sample,num_obs_carriers'")
    single_required.add_argument("-g", "--regional_freq_outfile", metavar='|',
                        help="File to store simulation output in csv format",
                        default=None)

    require_one = parser.add_argument_group('If climbing, require at least ' +\
                        'one of')
    require_one.add_argument("-H", "--hetfile", metavar='|',
                        help="File containing heterozygous probands")
    require_one.add_argument("-O", "--homfile", metavar='|',
                        help="File containing homozygous probands")

    batch_args = parser.add_argument_group('Args for batch validation.')
    batch_args.add_argument('-m', '--msub_validation',
                        help="Flag to perform batch validaiton on " +\
                        "a compute cluster using msub",
                        action="store_true")
    batch_args.add_argument('-F', '--calc_regional_freqs',
                        help="Flag to calculate expected regional allele " +\
                        "frequencies", action="store_true")
    batch_args.add_argument("-a", "--panel_file", metavar='|',
                        help="File containing panels for validation")
    batch_args.add_argument("-V", "--outpath", metavar='|',
                        help="Path to output validation results, each" +\
                        " panel in a timestamped directory")
    batch_args.add_argument("-R", "--resource_header", metavar='|',
                        help="File containing resource allocation header " +\
                        "to use when submitting jobs")
    batch_args.add_argument("-T", "--hom_het", metavar='|',
                        help="Indicate whether panels in --panel_file " +\
                        "are 'hom' or 'het' (not both)",
                        choices=['hom', 'het'])
    batch_args.add_argument("-t", "--test",
                        help="Write scripts but do not submit them",
                        action='store_true')

    args = parser.parse_args()

    if args.msub_validation is True:
        ## Check that proper validation args have been provided
        assert args.panel_file is not None, "Must provide validation panels"
        assert args.outpath is not None, "Must provide outpath"
        assert args.resource_header is not None, "Must provide resource " +\
                "allocation header when performing batch validation"

        ## TODO: Right now only one type of validation panel is supported +p4 id:164
        ## DONE: Check that this matches with constructed panels +p1 id:170

        import validation_msub
        validation_msub.main(args)
    else:
        assert args.climb_outfile is not None, "If not performing batch " +\
                "validation, climbing output file must be specified"
        assert args.control_outfile is not None, "If not performing batch " +\
                "validation, control likelihood output file must be specified"
        assert args.freq_params is not None, "If not performing " +\
                "batch validation, allele frequency parameters " +\
                "must be specified"
        main(args)
