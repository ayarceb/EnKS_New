#! /bin/bash

#============== Script for run the EnKS_MC for the LOTOS-EUROS  =======================
# Andres Yarce Botero-Santiago Lopez-Restrepo -Elias Niño
# Working from Barranquilla
#==========================================================================
#==================================================================================
#===Path where the program is running===
mydir='/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC/ENKS_MC_LOTOS'

#===Path LOTOS-EUROS MODEL (OJO carpeta de LEKF)===
LE='/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/lekf_4DEnVAR/lekf/v3.0.003-beta'

#===Path where netcdf-fortran and netcdf is installed===
OPT=${HOME}'/opt'
NETCDF_FORTRAN_HOME='/usr/lib64'
#NETCDF_HOME=${OPT}'/netcdf/4.4.0'
NETCDF_HOME=${OPT}'/home/dirac/miniconda3/pkgs/libnetcdf-4.7.3-nompi_h9f9fd6a_101'
#===Run ID====
runid='Prueba_numero_6_EnKS_MC'

if [ -f ${LE}/proj/eafit/000/rc/timerange.rc ]
then
	rm ${LE}/proj/eafit/000/rc/timerange.rc
	
fi	

#===Date of simulations====
if [ -d ${LE}/proj/eafit/000/rc/timerange.rc ]
then
	rm ${LE}/proj/eafit/000/rc/timerange.rc
fi	
start_date=20190202
echo 'timerange.start     :  2019-02-01 00:00:00'>>${LE}/proj/eafit/000/rc/timerange.rc
echo 'timerange.end       :  2019-02-05 00:00:00'>>${LE}/proj/eafit/000/rc/timerange.rc

#===Number of Ensembles===
Nens=40


#===Remove all temporal files====
if [ -d ${mydir}/temp ]
then
	rm ${mydir}/temp/*
fi
#===Path LOTOS-EUROS Ensemble Outputs===
LE_Outputs='/run/media/dirac/Datos/scratch/projects/'${runid}'/'${runid}'/output'

#==================================================================================
# 			      End Modified by user
#==================================================================================


#===Write runid, timerange and Ensemble files====
echo 'run.id             : '${runid}>${LE}/proj/eafit/000/rc/runid.rc
echo 'run.project             : '${runid}>${LE}/proj/eafit/000/rc/runproject.rc


echo 'kf.nmodes             : '${Nens}> ${LE}/proj/eafit/000/rc/N_Ense#mbles.rc
echo 'kf.nmodes             : '${Nens}> ${LE}/proj/eafit/000/rc/N_Ensembles.rc
echo ${LE_Outputs}>>${mydir}/temp/Ensembles.in
echo ${Nens}>>${mydir}/temp/Ensembles.in

#===Run Real State Simulation===
echo 'Running Model Real and Ensemble'

#====Run LOTOS-EUROS MODEL====  Comment if you dont want to run the ensemble
cd ${LE}

./launcher


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






