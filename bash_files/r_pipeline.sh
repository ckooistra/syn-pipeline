#!/bin/bash

# Get command line options and load flags into variables
# -n name that directory will be saved to if not provided subtype
# -s cancer subtype that we are looking at

while getopts ":n:s:t:d:" opt; do
	case ${opt} in
		n ) name=$OPTARG
			;;
		s ) subtype=$OPTARG
			;;
		d ) dir=$OPTARG
			;;
		\? ) echo "Invalid option ([-n] | [-s] | [-t])"
			;;
	esac
done

if [ -z $subtype ]; then
	echo "subtype [-s] required parameter"
	exit 1
fi

if [ -z $dir ]; then
	echo "output [-d] required parameter"
	exit 1
fi

if [ -z $name ]; then
	name=$subtype
fi

#R graphs folder path
r_graph_fp="/home/chris/Dropbox/BIN_3005/R_graphs/$dir/"

#Create Directory for analysis output if needed
if [ ! -d $r_graph_fp ]; then
	mkdir $r_graph_fp
fi	

data="/home/chris/Dropbox/BIN_3005/output.csv"
#Isolate data according to parameter entered
awk -F, -v subty=$subtype '($15==subty || NR==1 || $16==subty) {print $0}' $data > tmp

mkdir "$r_graph_fp$name"

Rscript '/home/chris/Dropbox/BIN_3005/Rscripts/graphs2.r' $subtype $name $dir 

rm tmp
