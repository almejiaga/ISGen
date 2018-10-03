#!/bin/bash
set -e

if [ "$#" -ne 7 ]; then
    echo "Bootstraps specified simulation results, first"
    echo "calculating symmetries within the pedigree,"
    echo "submitting each panel as a separate job via msub"
    echo
    echo "Usage: ./bootstrap_msub.py PEDFILE PANEL_FILES CLIMB_FILES"
    echo "CONTROL_FILES OUTPUT_FILES ITERATIONS RESOURCE_HEADER"
    exit 1
fi

pedfile=$1
mapfile -t panel_files < $2
mapfile -t climb_files < $3
mapfile -t control_files < $4
mapfile -t output_files < $5
iterations=$6
header_file=$7

## Check that each file had the same number of entries
if [ ${#panel_files[@]} != ${#climb_files[@]} ] || [ ${#panel_files[@]} != ${#control_files[@]} ] || [ ${#panel_files[@]} != ${#output_files[@]} ]; then
    echo "Input files not same length!"
    exit 1
fi

for ((i=0; i<${#panel_files[@]}; i++));	do
    ## If output file already exists, skip this iteration
    output_file=${output_files[$i]} 
    if test -e $output_file; then
        echo $output_file already exists - skipping
        continue
    fi

    ## Initialize script to submit
    panel_file=${panel_files[$i]} 
    script_file=$(dirname "${panel_file}")/bootstrap.sh
    cp $header_file $script_file
    echo '#PBS' -o $(dirname "${panel_file}")/ >> $script_file
    echo '#PBS' -e $(dirname "${panel_file}")/ >> $script_file
    echo >> $script_file
    echo 'cd $PBS_O_WORKDIR' >> $script_file
    echo cd ../bootstrap/ >> $script_file

    ## Calculate symmetries in the pedigree for the given panel inds
    symmetry_file=$(dirname "${panel_file}")/symmetry.txt
    echo python ped_symmetry.py -f $pedfile -p $panel_file -o $symmetry_file >> $script_file
    echo >> $script_file

    ## Perform bootstrapping
    climb_file=${climb_files[$i]} 
    control_file=${control_files[$i]} 
    echo python bootstrap.py single -l $climb_file -o $control_file -i $iterations -s $symmetry_file -c $output_file >> $script_file

    ## Submit job and pause so the scheduler doesn't get overloaded
    qsub $script_file
    sleep 1
done
