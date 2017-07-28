while read p; do
	/home/chris/Dropbox/BIN_3005/bash_files/r_pipeline.sh -s $p -n $p -d "histologies"
	echo $p
done < /home/chris/Dropbox/BIN_3005/python_files/supplemental/big_cancer_type.txt

