#!/bin/bash

# This is a script that runs all the samples and diffs them agains the ones in Samples/

RED='\033[1;31m'
NC='\033[0m' # No Color
compare () {
    cmp --silent ../Samples/$1 loc_$1 && echo "OK: ${1}"  || echo -e "${RED}Different${NC}: ${1}"
}

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
    time ../ezspectrum_code/ezFCF_linux.exe loc_$sample > loc_${sample}.out
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
    compare $spectrum
done

duschinsky="adenine.xml.spectrum_dushinsky cis_hcoh.xml.spectrum_dushinsky formaldehyde.xml.spectrum_dushinsky the_only_initial_state.xml.spectrum_dushinsky thymine.xml.spectrum_dushinsky"
for spectrum in $duschinsky
do
    compare $spectrum
done

echo ""
echo " = = = = = = = = = = "
echo " All tests completed "
echo " = = = = = = = = = = "
echo "" 

