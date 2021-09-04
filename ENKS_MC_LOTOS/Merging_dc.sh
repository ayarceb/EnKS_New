#! /bin/bash
runid='Prueba_numero_4_EnKS_MC'
#===Path LOTOS-EUROS Ensemble Outputs===
LE_Outputs='/run/media/dirac/Datos/scratch/projects/'${runid}'/'${runid}'/output'

start_date=20190201
#==Merging LE DC for each ensemble member ==
echo 'Merging LE Ensembles DC'
cd ${LE_Outputs}
let "j=0"
for i in $(ls LE_${runid}_dc_${start_date}_xi**a.nc)
	do
	let "j=j+1"
	if [ $j -lt 10 ]
	then
		ncks -O -h --mk_rec_dmn time LE_${runid}_dc_${start_date}_xi0${j}a.nc  Merge_x0${j}.nc
		mv LE_${runid}_dc_${start_date}_xi0${j}a.nc ..
		ncrcat -O -h Merge_x0${j}.nc LE_${runid}_dc_2*_xi0${j}a.nc Ens_dc_x0${j}.nc

	else
			ncks -O -h --mk_rec_dmn time LE_${runid}_dc_${start_date}_xi${j}a.nc  Merge_x${j}.nc
		mv LE_${runid}_dc_${start_date}_xi${j}a.nc ..
		ncrcat -O -h Merge_x${j}.nc LE_${runid}_dc_2*_xi${j}a.nc Ens_dc_x${j}.nc
	fi
	
	echo $j
	echo $i
	

done 

echo 'Merging Concentration surface files'
