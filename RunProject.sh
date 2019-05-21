#!/bin/bash
# Call function using bash
# Call function within folder with pictures.


for f in *.tif;
#do echo "Processing $f file.."; 
do if [[ $f == *"cy3"* ]]; then 
#	echo "Cy3 file fetched"
	f_2="${f/cy3/dapi}"
#	echo $f_2
#	echo "Rscript --vanilla /home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/Project.R --cy $f --dapi $f_2"

	Rscript --vanilla /home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/Project.R --cy $f --dapi $f_2 >> data.csv
	echo >> data.csv
#	echo "Finished, finding next file." 

else 
#	echo "DAPI file fetched, skipping"
	continue
fi

done
