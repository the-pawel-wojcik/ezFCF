#!/bin/bash

# This is a script cleans after the ./run.sh script

cleanProgress () {
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
        rm ${file}
    else
        echo ""
        echo "Cannot remove missing ${file}"
    fi
}


# Makse sure that dependencies are satisifed
if [ -e atomicMasses.xml ]
then 
    rm atomicMasses.xml
    echo " Local copy of atomicMasses.xml cleaned."
    echo ""
fi

# Clean a local copy of samples and print timing 
echo -n " Cleaning local samples and outputs"
samples="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml test.xml"
for sample in $samples
do
    cleanProgress $sample
    cleanProgress $sample "out"
done
echo " done."

echo -n " Cleaning parallel spectra"
parallel="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml test.xml"
for spectrum in $parallel
do
    cleanProgress $spectrum "spectrum_parallel"
done
echo " done."

echo -n " Cleaning duschinsky spectra"
duschinsky="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml"
for spectrum in $duschinsky
do
    cleanProgress $spectrum "spectrum_dushinsky"
done
echo " done."
