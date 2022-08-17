#!/bin/bash

# This is a script that runs all the samples and diffs them agains the ones in Samples/

# Makse sure that dependencies are satisifed
if [ ! -e atomicMasses.xml ]
then 
    cp ../atomicMasses.xml ./
    echo " Local copy of atomicMasses.xml created."
    echo ""
fi

echo ""
echo " = = = = = = = = = = = = = "
echo " Beginning the test "
echo " = = = = = = = = = = = = = "
echo "" 

# Run a local copy of samples and print timing 
samples="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml"
for sample in $samples
do
    echo ""
    echo $sample
    cp ../Samples/$sample ./loc_$sample
    time ../bin/ezFCF_linux.exe loc_$sample > loc_${sample}.out
done

echo ""
echo " = = = = = = = = = = = = = "
echo " Done running the examples "
echo " = = = = = = = = = = = = = "
echo "" 

# Diff parallel spectra
parallel="adenine.xml.spectrum_parallel cis_hcoh.xml.spectrum_parallel formaldehyde.xml.spectrum_parallel the_only_initial_state.xml.spectrum_parallel trans_hcoh_small.xml.spectrum_parallel trans_hcoh.xml.spectrum_parallel vg_phenolate.xml.spectrum_parallel"
for spectrum in $parallel
do
    # echo ${spectrum:0:-22}
    vimdiff loc_$spectrum ../Samples/$spectrum
done

duschinsky="adenine.xml.spectrum_dushinsky cis_hcoh.xml.spectrum_dushinsky formaldehyde.xml.spectrum_dushinsky the_only_initial_state.xml.spectrum_dushinsky thymine.xml.spectrum_dushinsky"
for spectrum in $duschinsky
do
    # echo ${spectrum:0:-22}
    vimdiff loc_$spectrum ../Samples/$spectrum
done
