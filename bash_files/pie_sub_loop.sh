while read p; do
	Rscript /home/chris/Dropbox/BIN_3005/code/Rscripts/pie_chart_outline.r $p 
	echo $p
done < /home/chris/Dropbox/BIN_3005/code/supplemental/cancer_spots.txt


