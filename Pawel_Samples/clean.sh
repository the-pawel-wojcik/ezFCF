#!/bin/bash

# This is a script cleans after the ./run.sh script

# Makse sure that dependencies are satisifed
if [ -e atomicMasses.xml ]
then 
    rm atomicMasses.xml
    echo " Local copy of atomicMasses.xml cleaned."
    echo ""
fi

# Clean a local copy of samples and print timing 
echo -n " Cleaning local samples and outputs"
samples="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml"
for sample in $samples
do
    echo -n "."
    if [ -e $sample ]
    then
        rm $sample
    fi
    if [ -e ${sample}.out ]
    then
        rm ${sample}.out
    fi
done
echo " done."

echo -n " Cleaning parallel spectra"
parallel="adenine.xml.spectrum_parallel cis_hcoh.xml.spectrum_parallel formaldehyde.xml.spectrum_parallel the_only_initial_state.xml.spectrum_parallel trans_hcoh_small.xml.spectrum_parallel trans_hcoh.xml.spectrum_parallel vg_phenolate.xml.spectrum_parallel"
for spectrum in $parallel
do
    echo -n "."
    if [ -e $spectrum ]
    then
        rm $spectrum
    fi
done
echo " done."

echo -n " Cleaning duschinsky spectra"
duschinsky="adenine.xml.spectrum_dushinsky cis_hcoh.xml.spectrum_dushinsky formaldehyde.xml.spectrum_dushinsky the_only_initial_state.xml.spectrum_dushinsky thymine.xml.spectrum_dushinsky"
for spectrum in $duschinsky
do
    echo -n "."
    if [ -e $spectrum ]
    then
        rm $spectrum
    fi
done
echo " done."
