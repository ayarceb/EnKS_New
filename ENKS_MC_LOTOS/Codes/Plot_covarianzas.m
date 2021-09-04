
% code to plot the covariances from the B cholesky matrices

clc;clear all;
close all
%load cholesky_matrices.mat
load B_state_dc_lag1.mat


mydir='/run/media/dirac/Datos/scratch/projects/Prueba_numero_4_EnKS_MC/Prueba_numero_4_EnKS_MC/output';
cd(mydir)
%  subplot(2,3,1);imagesc(B_state_dc_lag0);colorbar;caxis([0 1e-15]);colormap(jet);title('L 0');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
%  subplot(2,3,2);imagesc(B_state_dc_lag1);colorbar;caxis([0 1e-22]);colormap(jet);title('L 1');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
%  subplot(2,3,3);imagesc(B_state_dc_lag2);colorbar;caxis([0 1e-23]);colormap(jet);title('L 2');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
%  subplot(2,3,4);imagesc(B_state_dc_lag3);colorbar;caxis([0 1e-24]);colormap(jet);title('L 3');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
%  subplot(2,3,5);imagesc(B_state_dc_lag4);colorbar;caxis([0 1e-25]);colormap(jet);title('L 4');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')

lat=ncread('LE_Prueba_numero_4_EnKS_MC_column_20190205_xi02a.nc','latitude');
lon=ncread('LE_Prueba_numero_4_EnKS_MC_column_20190205_xi02a.nc','longitude');

mydir='/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC/ENKS_MC_LOTOS/Codes';
cd(mydir)



figure
for i=1:100
Sample=B_state_dc_lag1(2800+i,:);

Sample_=normalize(reshape(Sample(2342:end-1),[45,52]));   %  part to reorganize the cholesky matrices


imagesc(lon,lat,Sample_'); colormap jet;set(gca,'YDir','normal'); hold on

S=shaperead('/home/dirac/Dropbox/2020/4DENVAR_PAPER/Prueba4DEnVAR_25ensembles/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.9]);ylim([9 13.5]);h=colorbar;ylabel(h,'1e15 mlc/cm2');xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')


colorbar;caxis([0 1])
pause(0.1)
end
