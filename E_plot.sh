#!/bin/bash -f

filename=$1		#initial structure xyz file
bond=1.31 		#bond length of initial structure

listxscale=" "
for ((xscale=950;xscale<=1050;xscale++))
do
	xscale_dec=$( echo "scale=3;($xscale/1000)" | bc)
	listxscale="$listxscale $xscale_dec"
done

listyscale=" "
for ((yscale=950;yscale<=1050;yscale++))
do
	yscale_dec=$( echo "scale=3;($yscale/1000)" | bc)
	listyscale="$listyscale $yscale_dec"
done

listzscale=" "
for ((zscale=950;zscale<=1050;zscale++))
do
	zscale_dec=$( echo "scale=3;($zscale/1000)" | bc)
	listzscale="$listzscale $zscale_dec"
done

#Calculate energies at each bond*xscale
for xscale in $listxscale
do
#Scale input structure
cd xyzfiles
./scalecell ../$filename scaled.xyz $xscale 1 1
./genkgrid scaled.xyz temp.kpts 10 1 1
cd ..

./singleE_main xyzfiles/scaled.xyz xyzfiles/temp.kpts > energy
grep "Etot" energy > temp
awk '{printf "%.6f\n", (($3))}' temp >> energies
echo $( echo "scale=6;($xscale*$bond)" | bc ) >> bonds
done

#Combine bond lengths and energies into final output file
paste -d"\t" bonds energies > E_vs_bond.dat

#Clean up
rm energy temp energies bonds 

