#!/bin/bash

age=0
	while (( $age <= 140 ))
	do

		step=`echo "scale=1; $age + 1" | bc`
	#mv -v [${age}Ma_Golonka.xyz]_Hyd_zonal_185_${step}.vtk [${age}Ma_Golonka.xyz]_Atm_zonal.vtk
	
	mv -v [${age}Ma_Golonka.xyz]_Atm_longal_45_${step}.vtk [${age}Ma_Golonka.xyz]_Atm_longal.vtk
		
age=$(($age + 1))
done