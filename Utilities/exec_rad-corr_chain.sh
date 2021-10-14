#!/bin/bash

#######################################################
#                                                     #
#  Script made to execute the chain of programs       #
#        to calculate electron radiative corrections  #
#                                                     #
#######################################################

function print_help() {
    echo "SCRIPT: exec_rad-corr_chain.sh"
    echo "======================="
    echo "./exec_rad-corr_chain.sh --part <part>"
    echo "Where:"
    echo "  <part>  = selects particle (eta, omega)"
    echo "Example:"
    echo "  ./exec_rad-corr_chain.sh --part omega"
    exit
}

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--part" ]]; then
            part=${arr[$((ic+1))]}
        else
            echo "ERROR: unrecognized argument: ${arr[$((ic))]}.";
            print_help;
        fi
        ((ic+=2))
    done
}

function print_args() {
    echo "SCRIPT: exec_rad-corr_chain.sh"
    echo "=============================="
    echo "part = ${part}"
}

################
###   Main   ###
################

if [[ -z "${EXTERNALS}" ]]; then
    echo "ERROR: variable EXTERNALS is unset."
    exit 1
fi

if [[ -z "${PRODIR}" ]]; then
    echo "ERROR: variable PRODIR is unset."
    exit 1
fi

if [[ ${#} -ne 2 ]]; then
    echo "ERROR: ${#} arguments were provided, they should be 2."
    print_help
fi

argArray=("$@")
process_args "${argArray[@]}"
print_args

targets=("D" "C" "Fe" "Pb")
targets_v2=("d2" "C12" "Fe56" "Pb208")

# move to Utilities/
cd ${EXTERNALS}/Utilities
make clean; make
cd ${EXTERNALS}/Utilities/bin
./GetBinning -p${part}

# loop over targets
for ((t=0; t<4; t++)); do
    # execute GetCentroids
    ./GetCentroids -t${targets[$t]} -p${part}

    # copy centroids info to RUNPLAN
    cp -v centroids_${part}_${targets[$t]}.txt ${EXTERNALS}/RUNPLAN/

    # enter RUNPLAN, modify make_plan.py to open desired centroids file
    cd ${EXTERNALS}/RUNPLAN/
    sed -i "s|f=open(\"centroids.txt\");|f=open(\"centroids_${part}_${targets[$t]}.txt\");|g" make_plan.py

    # execute it
    python make_plan.py > clas_kin.inp

    # move to main dir, compile EXTERNALS and execute it
    cd ${EXTERNALS}
    make clean; make
    ./run_extern.sh clas${targets_v2[$t]}

    # restore make_plan.py
    cd ${EXTERNALS}/RUNPLAN/
    sed -i "s|f=open(\"centroids_${part}_${targets[$t]}.txt\");|f=open(\"centroids.txt\");|g" make_plan.py

    # go back to Utilities/bin
    cd ${EXTERNALS}/Utilities/bin
done
