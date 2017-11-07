#!/bin/bash


while getopts ":s:t:" opt; do
	case ${opt} in
		s ) subtype=$OPTARG
			;;
		t ) codon=$OPTARG
			;;
		\? ) echo "Invalid option ([-d] | [-s])"
			;;
	esac
done

fp="/home/chris/Dropbox/BIN_3005/mutation_info/$subtype"

if [ ! -d $fp ]; then
	mkdir $fp
fi

cd $fp
/home/chris/anaconda3/bin/python3 /home/chris/Dropbox/BIN_3005/code/python/codons.py -s $subtype -t $codon 
