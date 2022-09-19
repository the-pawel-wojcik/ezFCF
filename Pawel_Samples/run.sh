#!/bin/bash

# This script runs locally copies of all inputs from Samples/ and diffs 
# the results against the original outputs.

# Specify path and filenme of ezFCF executable
ezFCF="../ezspectrum_code/ezFCF_linux.exe" # use this one if you compiled the program
# ezFCF="../bin/ezFCF_linux.exe"

# Check if ezFCF file exists and is executable
if [ ! -x $ezFCF ]
then
    echo "Error. The ezFCF exacutable at $ezFCF does not exist."
    echo "Make sure the file exists and is executable (chmod +x)"
    echo " or fix the path and filename inside this scirpt."
    exit 1
fi

# Make sure that dependencies are satisifed
if [ ! -e atomicMasses.xml ]
then 
    cp ../atomicMasses.xml ./
    echo " Local copy of atomicMasses.xml created."
    echo ""
fi

echo ""
echo " = = = = = = = = = = = = = = ="
echo "  Running sample calculations "
echo " = = = = = = = = = = = = = = ="
echo "" 

# Run a local copy of samples and print timing 
samples="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml"
for sample in $samples
do
    echo ""
    echo $sample
    cp ../Samples/$sample ./$sample
    time $ezFCF $sample > ${sample}.out
done

echo ""
echo " = = = = = = = = = ="
echo "  Samples are ready"
echo " = = = = = = = = = ="
echo "" 

RED='\033[1;31m'
NC='\033[0m' # No Color
compare () {
    cmp --silent ../Samples/$1 $1 && echo "OK: ${1}"  || echo -e "${RED}Different${NC}: ${1}"
}

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
echo " = = = = = = = = = = ="
echo "  Comparison complete "
echo " = = = = = = = = = = ="
echo "" 
