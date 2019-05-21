#!/bin/bash
# Call function using bash
# Call function within folder with pictures.#Run with ls | parallel bash RunProject_Parallel.sh
# ls | parallel bash /home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/RunProject_parallel.sh


if [[ $1 == *"cy3"* ]]; then
	echo "Cy3 file fetched"
	f_2="${1/cy3/dapi}"
	echo $f_2
#	echo "Rscript --vanilla /home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/Project.R --cy $f --dapi $f_2"
	Rscript --vanilla /home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/Project.R --cy $1 --dapi $f_2 >> data_parallel.csv
	echo >> data_parallel.csv
#	echo "Finished, finding next file." 
else
	echo "Not a Cy3 file"
fi




#doit() {
#    f=*
#    echo $f
#    f_2="${f/cy3/dapi}"
#    Rscript --vanilla /home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/Project.R --cy $f --dapi $f_2 >> data_parallel.csv
#    echo >> data_parallel.csv
#}
#export -f doit
#parallel -j 15 doit ::: *cy3*.tif




