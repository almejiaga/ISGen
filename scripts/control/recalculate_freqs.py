import sys, os
sys.path.append('/Users/dnelson/project/anc_finder/scripts/control')
import numpy as np
import tables
import regional_freq as rf

try:
    freq_params = sys.argv[1]
    climb_file = sys.argv[2]
    control_file = sys.argv[3]
    merged_file = sys.argv[4]
except:
    print "Usage: python recalculate_freqs.py POP_SIZE,SAMPLE_SIZE,NUM_OBS_ALLELES CONTROL_FILE CLIMB_FILE FREQ_FILE"
    raise

with tables.open_file(merged_file, 'r') as f:
    nodes = f.walk_nodes(f.root)
    rootnode = next(nodes)
    regions = [node.name for node in nodes if node.name != 'regional_expected_freqs']

with tables.open_file(merged_file, 'r') as f:
        node = f.get_node(f.root, 'All Probands')
        global_expectations = node[:]
        ntrees = global_expectations.shape[0]
        
with tables.open_file(climb_file, 'r') as f:
    climb_liks = f.root.liks[:ntrees]
    
with tables.open_file(control_file, 'r') as f:
    control_liks = f.root.control_liks[:ntrees]

    
regional_expectations = []

for region in regions:
    # print "Calculating expected allele frequencies for region", region

    with tables.open_file(merged_file, 'r') as f:
        node = f.get_node(f.root, region)
        regional_freqs = node[:]
        
    numerator = 0.
    denominator = 0.

    for E_r, c, t, E_F in zip(regional_freqs, control_liks, climb_liks,
                            global_expectations):
        N, n, k = [int(x.strip()) for x in freq_params.split(',')]
        F_obs = 1. * k / n * N

        E_rF = rf.get_E_rF(E_F, E_r, freq_params)
        numerator += E_rF * c * t
        denominator += c * t
        
    if numerator == 0:
        ## Avoids division by zero if expected frequency is zero
        regional_expectations.append(0)
    else:
        assert denominator != 0
        regional_expectation = numerator / denominator
        regional_expectations.append(regional_expectation)

for region, expectation in zip(regions, regional_expectations):
    print region + ',' + str(expectation)

