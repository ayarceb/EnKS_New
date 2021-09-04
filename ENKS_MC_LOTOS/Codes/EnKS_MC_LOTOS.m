% Code to implement the EnKS_MC technique with the LOTOS-EUROS

clc;close all;clear all

% Variable characteristics    Size no2:       48x50x1x97
%                             Dimensions: longitude,latitude,level,time
%                             Datatype:   single
%                    Attributes:
%                       standard_name   = 'mole_fraction_of_nitrogen_dioxide_in_air'
%                       long_name       = 'volume mixing ratio of NO2 in humid air'
%                       units           = 'mole mole-1'
%                       _CoordinateAxes = 'time level latitude longitude'
%                       molemass        = 0.046005
%                       molemass_unit   = 'kg mole-1'


% Coordenadas        (latitude,longitude)
% Barranquilla       Latitud: 10.9878,    Longitud: -74.7889      (22,11)
% Santa Marta        Latitud: 11.24079,   Longitud: -74.19904      (25,17)
% Cartagena          Latitud: 10.3884731, Longitud: -75.488593     (16,3)
% Mina               Latitud: 9.5829857,  Longitud: -73.479376      (7,25)
% Valledupar         Latitud: 10.4646054, Longitud: -73.2599952      (17,28)

ciudades={'Barranquilla ','Santa Marta ','Cartagena ','Mina Drummond ','Valledupar '};

lon=[11,17,3,25,28];lat=[22,25,16,7,17];

name_run='Prueba_numero_4_EnKS_MC';

%%===Path where the output are located from the EnKS_MC ensemble first propagation
mydir=append('/run/media/dirac/Datos/scratch/projects/',name_run,'/',name_run,'/output');cd(mydir)

system('rm no2_column_ens*.nc');system('rm Merge_*.nc'); system('rm Ens_dc_*.nc');
cd ..
system('mv LE_Prueba_numero_4_EnKS_MC_dc_*.nc output')
mydir=append('/run/media/dirac/Datos/scratch/projects/',name_run,'/',name_run,'/output');cd(mydir)


no2=zeros(45,52,97,40); no2_dc=zeros(45,52,97,40);

t1 = datetime(2019,2,1,0,0,0);t2 = datetime(2019,2,5,0,0,0);t = t1:hours(1):t2;


ens=40; %NO2_sfc Ensemble number concatenate
for j=1:5
for  n=1:40  %cycle to concatenate the LOTOS-EUROS output, concatenate and read 
if n<10    
system(append(sprintf(append('cdo select,name=no2  LE_',name_run,'_conc-sfc_2019020*_xi0%ia.nc no2'),n),sprintf('_column_ens_%i.nc',n)))  % Concatenate all ensemble member 
no2(:,:,:,n)=squeeze(ncread(sprintf('no2_column_ens_%i.nc',n),'no2'));
NO2_state(:,:,:,n)=no2(:,:,:,n);
end
 
if n>=10    
system(append(sprintf(append('cdo select,name=no2 LE_',name_run,'_conc-sfc_2019020*_xi%ia.nc no2'),n),sprintf('_column_ens_%i.nc',n)))  % Concatenate all ensemble member 
no2(:,:,:,n)=squeeze(ncread(sprintf('no2_column_ens_%i.nc',n),'no2'));
NO2_state(:,:,:,n)=no2(:,:,:,n);
end

latitude=ncread('no2_column_ens_1.nc','latitude'); 
longitude=ncread('no2_column_ens_1.nc','longitude');

no2_=squeeze(no2(lon(j),lat(j),5:end,n));
subplot(2,3,j)
h=plot(t(1:end-4),no2_,'-b','LineWidth',2);
hold on

mycolor = [224/256 224/256 224/256];
set(h,'Color', mycolor);
end

no2_mean=squeeze(mean(squeeze(no2(lon(j),lat(j),:,:))')); 
p=plot(t(1:end-4),no2_mean(5:end),'k','LineWidth',3);
grid on
ylabel(sprintf('Mole fraction of nitrogen dioxide \n in air mole mole-1'));
title(append(ciudades{j}, sprintf('lat= %1.2f ',latitude(j)),'°', sprintf('lon= %1.2f ',longitude(j)),'°'));

legend([h,p],'Ensemble members','Mean');
end

mydir=append('/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC/ENKS_MC_LOTOS');cd(mydir)

system('./Merging_dc.sh')  % Bash code to merge dc factors

mydir=append('/run/media/dirac/Datos/scratch/projects/',name_run,'/',name_run,'/output');cd(mydir)


%% PLOT DC

figure

for j=1:5
for  n=1:ens  %cycle to concatenate the LOTOS-EUROS output, concatenate and read 
if n<10    
no2_dc(:,:,:,n)=squeeze(ncread(sprintf('Ens_dc_x0%i.nc',n),'dc'));
end
 
if n>=10    
no2_dc(:,:,:,n)=squeeze(ncread(sprintf('Ens_dc_x%i.nc',n),'dc'));
end

latitude=ncread('no2_column_ens_1.nc','latitude'); 
longitude=ncread('no2_column_ens_1.nc','longitude');

no2_dc_=squeeze(no2_dc(lon(j),lat(j),5:end,n));
subplot(2,3,j)
h=plot(t(1:end-4),no2_dc_,'-b','LineWidth',2);
hold on

mycolor = [224/256 224/256 224/256];
set(h,'Color', mycolor);
end

no2_mean_dc=squeeze(mean(squeeze(no2_dc(lon(j),lat(j),:,:))')); 
p=plot(t(1:end-4),no2_mean_dc(5:end),'k','LineWidth',3);
grid on
ylabel(sprintf('DC factor'));
title(append(ciudades{j}, sprintf('lat= %1.2f ',latitude(j)),'°', sprintf('lon= %1.2f ',longitude(j)),'°'));


end

figure(1)
Xb=nan(45,52,97,40);
Fe=nan(45,52,97,40);

Xb_mean=nan(5,97);
Fe_mean=nan(5,97);

% for j=1:5   % geolocalizaciones
% for l=1:4   % lags
% Xb(:,:,18+(-1+l)*24,:)=no2(:,:,18+(-1+l)*24,:);  % Background  states
% Fe(:,:,18+(-1+l)*24,:)=no2_dc(:,:,18+(-1+l)*24,:);   % Background parameters
% Xb_mean(j,18+(-1+l)*24)=squeeze(mean(squeeze(no2(lon(j),lat(j),18+(-1+l)*24,:))'));  % Background mean states
% Fe_mean(j,18+(-1+l)*24)=squeeze(mean(squeeze(no2_dc(lon(j),lat(j),18+(-1+l)*24,:))'));   % Background mean parameters
% end
% end
%     
% 
% for j=1:5
% subplot(2,3,j)
%     hh=plot(t(1:end-4),squeeze(Xb(lon(j),lat(j),5:end,1)),'Marker', 'o', 'MarkerSize', 4, ...
%      'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     pp=plot(t(1:end-4),squeeze(Xb_mean(j,5:end)),'Marker', 'o', 'MarkerSize', 4, ...
%      'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% legend([h,p,hh,pp],'Ensemble members','Mean','Background ensemble member state', 'Background mean state')
% end
% 
% 
% 
% figure(2)
% for j=1:5
% subplot(2,3,j)
%    hhh= plot(t(1:end-4),squeeze(Fe(lon(j),lat(j),5:end,1)),'Marker', 'o', 'MarkerSize', 4, ...
%       'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%   ppp= plot(t(1:end-4),squeeze(Fe_mean(j,5:end)),'Marker', 'o', 'MarkerSize', 4, ...
%      'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%  legend([h,p,hhh,ppp],'Ensemble members','Mean','Background ensemble member state', 'Background mean state')
% end