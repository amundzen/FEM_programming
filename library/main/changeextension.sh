for FN in *.f95
do
	mv "${FN}" "${FN%f95}f90"
done