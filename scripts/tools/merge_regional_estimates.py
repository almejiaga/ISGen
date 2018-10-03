import tables
import numpy as np
import sys, os

if len(sys.argv) != 3:
    print("""Usage:
    
    python merge_output.py regional_freq_file_paths merged_freq_file
    """)
    sys.exit()

freq_file_paths = sys.argv[1]
freq_merge_file = sys.argv[2]

freq_files = sorted([line.strip() for line in open(freq_file_paths)])

with tables.open_file(freq_files[0], 'r') as f:
    nodes = [n for n in f.walk_nodes() if n._v_name != '/' and n._v_name != 'regional_expected_freqs']
    region_names = [n._v_name for n in nodes]
    region_sizes = [n._v_attrs['ninds'] for n in nodes]

with tables.open_file(freq_merge_file, 'w') as f:
    node_dict = {}
    for name, size in zip(region_names, region_sizes):
        node_dict[name] = tables.EArray(f.root,
            name, tables.FloatAtom(), shape=(0,))
        node_dict[name]._v_attrs['ninds'] = size

    for fname in freq_files:
        print("Merging", fname)
        with tables.open_file(fname, 'r') as g:
            for name in region_names:
                print("Merging", name)
                new_node = g.get_node(g.root, name)
                node_dict[name].append(new_node[:])

# with tables.open_file(climb_merge_file, 'w') as climb_merge_file:
#     ancs = tables.EArray(climb_merge_file.root, 'ancs', tables.IntAtom(), shape=(0,))
#     genotypes = tables.VLArray(climb_merge_file.root, 'genotypes', tables.IntAtom())
#     liks = tables.EArray(climb_merge_file.root, 'liks', tables.IntAtom(), shape=(0,))
#     trees = tables.VLArray(climb_merge_file.root, 'trees', tables.IntAtom())
#
#     for fname in climb_files:
#         print("Merging", fname)
#         with tables.open_file(fname, 'r') as f:
#             ancs.append(f.root.ancs[:])
#             liks.append(f.root.liks[:])
#
#             for line in f.root.trees[:]:
#                 trees.append(line)
#
#             for line in f.root.genotypes[:]:
#                 genotypes.append(line)
