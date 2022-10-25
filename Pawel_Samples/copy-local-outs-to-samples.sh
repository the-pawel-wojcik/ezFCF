#!/bin/bash

# This is a script cleans after the ./run.sh script

Progress () {
    # Prepare filename
    file=$1
    if [ ! -z $2 ]
    then
        file=${1}.${2}
    fi
    # remove if the file exists and print progress bar
    if [ -e ${file} ]
    then
        echo -n "."
        cp ${file} ../Samples/
    else
        echo ""
        echo "Cannot copy missing ${file}"
    fi
}

# Copy the local output of samples to ../Samples/
echo -n " Overwriting Samples with local outputs"
samples="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml test.xml"
for sample in $samples
do
    Progress $sample "out"
done
echo " done."

echo -n " Overwriting Samples with local parallel spectra"
parallel="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml test.xml"
for spectrum in $parallel
do
    Progress $spectrum "spectrum_parallel"
done
echo " done."

echo -n " Overwriting Samples with local duschinsky spectra"
duschinsky="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml"
for spectrum in $duschinsky
do
    Progress $spectrum "spectrum_dushinsky"
done
echo " done."
