import sys, os
import numpy as np
import pandas as pd
import tables
from collections import defaultdict
sys.path.append('/Users/dnelson/project/anc_finder/scripts/control')
sys.path.append('/Users/dnelson/project/anc_finder/scripts/bootstrap')
import bootstrap
import regional_freq as rf


try:
    freq_params = sys.argv[1]
    climb_file = sys.argv[2]
    control_file = sys.argv[3]
    merged_file = sys.argv[4]
    outfile = os.path.expanduser(sys.argv[5])
    iterations = int(sys.argv[6])
except:
    print "Usage: python recalculate_freqs.py POP_SIZE,SAMPLE_SIZE,NUM_OBS_ALLELES CLIMB_FILE CONTROL_FILE FREQ_FILE OUTPUT_FILE" + \
            " BOOTSTRAP_ITER"
    print ""
    raise

print "Performing", iterations, "bootstrap resamples"

N, n, k = [int(x.strip()) for x in freq_params.split(',')]
F_obs = 1. * k / n * N

data = dict()

with tables.open_file(merged_file, 'r') as f:
    nodes = f.walk_nodes(f.root)
    rootnode = next(nodes)
    regions = [node.name for node in nodes if node.name != 'regional_expected_freqs']

with tables.open_file(merged_file, 'r') as f:
    node = f.get_node(f.root, 'All Probands')
    data['global'] = node[:]
        
with tables.open_file(climb_file, 'r') as f:
    data['climb_liks'] = f.root.liks[:]
    
with tables.open_file(control_file, 'r') as f:
    data['control_liks'] = f.root.control_liks[:]

ninds_dict = {}
with tables.open_file(merged_file, 'r') as f:
    for region in regions:
        print "Loading region", region
        data[region] = []
        node = f.get_node(f.root, region)
        ninds_dict[region] = node._v_attrs['ninds']
        data[region] = node[:]

df = pd.DataFrame.from_dict(data)
df['total_liks'] = df[['climb_liks', 'control_liks']].sum(axis=1).map(lambda x: 2.**x)
total_lik_sum = df['total_liks'].sum()
bootresults = []

for i in xrange(iterations):
    if i % 100 == 0:
        print "Iteration", i, "of", iterations
    resamp = df.sample(frac=1, replace=True)
    result = resamp[regions].multiply(resamp['total_liks'].values, axis=0)
    result = result.sum() / total_lik_sum

    ## Approximation to condition on global allele frequency
    correction = F_obs / result['All Probands']
    result = result * correction

    ## Store allele frequency
    for region in regions:
        result[region] = result[region] / ninds_dict[region]

    result.reset_index(level=0)
    bootresults.append(result)

bootresults = pd.DataFrame.from_dict(map(dict,bootresults))
bootstrap.Bootstrap.save_boot(bootresults, outfile)
