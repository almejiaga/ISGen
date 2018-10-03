import tables
import numpy as np
import sys, os

if len(sys.argv) != 5:
    print("""Usage:
    
    python merge_output.py climb_file_paths merged_climb_file control_file_paths merged_control_file
    """)
    sys.exit()

climb_file_paths = sys.argv[1]
climb_merge_file = sys.argv[2]
control_file_paths = sys.argv[3]
control_merge_file = sys.argv[4]

climb_files = sorted([line.strip() for line in open(climb_file_paths)])
control_files = sorted([line.strip() for line in open(control_file_paths)])

assert len(climb_files) == len(control_files)

with tables.open_file(control_merge_file, 'w') as control_merge_file:
    control_liks = tables.EArray(control_merge_file.root,
            'control_liks', tables.IntAtom(), shape=(0,))
    tot_liks = tables.EArray(control_merge_file.root,
            'tot_liks', tables.IntAtom(), shape=(0,))

    for fname in control_files:
        print("Merging", fname)
        with tables.open_file(fname, 'r') as f:
            control_liks.append(f.root.control_liks[:])
            tot_liks.append(f.root.tot_liks[:])

with tables.open_file(climb_merge_file, 'w') as climb_merge_file:
    ancs = tables.EArray(climb_merge_file.root, 'ancs', tables.IntAtom(), shape=(0,))
    genotypes = tables.VLArray(climb_merge_file.root, 'genotypes', tables.IntAtom())
    liks = tables.EArray(climb_merge_file.root, 'liks', tables.IntAtom(), shape=(0,))
    trees = tables.VLArray(climb_merge_file.root, 'trees', tables.IntAtom())

    for fname in climb_files:
        print("Merging", fname)
        with tables.open_file(fname, 'r') as f:
            ancs.append(f.root.ancs[:])
            liks.append(f.root.liks[:])

            for line in f.root.trees[:]:
                trees.append(line)

            for line in f.root.genotypes[:]:
                genotypes.append(line)
