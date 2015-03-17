#!/bin/bash -f

filename=$1		#initial structure xyz file
bond=1.42 		#bond length of initial structure

listxscale=" "
for ((xscale=950;xscale<=1050;xscale++))
do
	xscale_dec=$( echo "scale=3;($xscale/1000)" | bc)
	listxscale="$listxscale $xscale_dec"
done

#Calculate energies at each bond*xscale
for xscale in $listxscale
do
#Scale input structure
cd xyzfiles
./scalecell ../$filename scaled.xyz $xscale $xscale 1
./genkgrid scaled.xyz temp.kpts 3 3 1
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

