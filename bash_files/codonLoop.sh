#!/bin/bash


while getopts ":t:" opt; do
	case ${opt} in
		t ) codon=$OPTARG
			;;
		\? ) echo "Invalid option ([-d] | [-s])"
			;;
	esac
done
while read p; do
	/home/chris/Dropbox/BIN_3005/code/bash_files/codon.sh -s $p -t $codon
	echo $p
done < /home/chris/Dropbox/BIN_3005/code/supplemental/cancer_places.txt
