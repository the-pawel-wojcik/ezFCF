#!/bin/bash

# This script runs locally copies of all inputs from Samples/ and diffs 
# the results against the original outputs.

# Specify path and filenme of ezFCF executable
ezFCF="../ezspectrum_code/ezFCF_linux.exe" # use this one if you compiled the program
# ezFCF="../bin/ezFCF_linux.exe"

if [[ $1 == "clean" ]]
then
    ./clean.sh
    cd ../ezspectrum_code
    make clean
    make
    cd ../Pawel_Samples
fi

if [[ $1 == "make" ]]
then
    ./clean.sh
    cd ../ezspectrum_code
    make
    cd ../Pawel_Samples
fi

# Check if ezFCF file exists and is executable
if [ ! -x $ezFCF ]
then
    echo "Error."
    echo "  Missing executable ezFCF at $ezFCF."
    echo "  Make sure the file exists and is executable (chmod +x)"
    echo "  or fix the path and filename inside this scirpt."
    exit 1
fi

# Make sure that dependencies are satisifed
if [ ! -e atomicMasses.xml ]
then 
    cp ../atomicMasses.xml ./
    echo -e " Local copy of atomicMasses.xml created.\n"
fi

echo ""
echo " = = = = = = = = = = = = = = ="
echo "  Running sample calculations "
echo " = = = = = = = = = = = = = = ="
echo "" 

# Run a local copy of samples and print timing 
samples="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml"
# samples=""
for sample in $samples
do
    echo ""
    echo $sample
    cp ../Samples/$sample ./$sample
    time $ezFCF $sample > ${sample}.out
done

# run an extra test from InputScripts
echo ""
echo test.xml
cp ../InputScripts/test.xml ./
time $ezFCF test.xml > test.xml.out

echo ""
echo " = = = = = = = = = ="
echo "  Samples are ready"
echo " = = = = = = = = = ="
echo "" 

RED='\033[1;31m'
NC='\033[0m' # No Color
compare () {
    directory=Samples
    if [ ! -z $2 ]
    then
        directory=$2
    fi
    cmp --silent ../${directory}/$1 $1 && echo "OK: ${1}"  || echo -e "${RED}Different${NC}: ${1}"
}

# Diff parallel spectra
parallel="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml trans_hcoh_small.xml trans_hcoh.xml vg_phenolate.xml"
for spectrum in $parallel
do
    compare ${spectrum}.spectrum_parallel
done

duschinsky="adenine.xml cis_hcoh.xml formaldehyde.xml the_only_initial_state.xml thymine.xml"
for spectrum in $duschinsky
do
    compare ${spectrum}.spectrum_dushinsky
done

compare test.xml.spectrum_parallel InputScripts

echo ""
echo " = = = = = = = = = = ="
echo "  Comparison complete "
echo " = = = = = = = = = = ="
echo "" 
