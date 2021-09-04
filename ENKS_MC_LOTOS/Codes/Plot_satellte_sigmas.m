% code to plot the satellite files from the TROPOMI simulation of the
% EnKS-MC paper

clc;close all;clear all

mydir='/run/media/dirac/Datos/scratch/projects/Prueba_numero_4_EnKS_MC/Prueba_numero_4_EnKS_MC/output';
cd(mydir)



date={'2019-02-01','2019-02-02','2019-02-03','2019-02-03','2019-02-04'};
tit={'NO_2 LOTOS-EUROS Column','Vertical column density','Simulated retrieval'};

files={'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190201_1800.nc','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190202_1800.nc',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1800.nc','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1900.nc',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190204_1900.nc'};

filesbase={'LE_Prueba_numero_4_EnKS_MC_column_20190201_xb.nc','LE_Prueba_numero_4_EnKS_MC_column_20190202_xb.nc',...
    'LE_Prueba_numero_4_EnKS_MC_column_20190203_xb.nc','LE_Prueba_numero_4_EnKS_MC_column_20190203_xb.nc',...
    'LE_Prueba_numero_4_EnKS_MC_column_20190204_xb.nc'};

random_n=[1000,1500,2000];

%  ilat=ncread(files{1},'ilat')
%  n=length(ilat)
%  randam_obs_1=randi([1 n],1,random_n(1))
%  ilat=ncread(files{2},'ilat')
%  n=length(ilat)
%  randam_obs_2=randi([1 n],1,random_n(2))
%  ilat=ncread(files{3},'ilat')
%  n=length(ilat)
%  randam_obs_3=randi([1 n],1,random_n(3))
%  ilat=ncread(files{4},'ilat')
%  n=length(ilat)
%  randam_obs_4=randi([1 n],1,random_n(4))
%  ilat=ncread(files{5},'ilat')
%  n=length(ilat)
%  randam_obs_5=randi([1 n],1,random_n(3))

cont=1;
for k=1:3
for i=1:5
    
if ((cont~=4)&&(cont~=9)&&(cont~=14))
%     figure(i)
i
k
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
   
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);
    
    randam_obs{k,i,:}=randi([1 n],1,random_n(k));
  
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{k,i,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
        end
   end 
    end
if (cont<=4)   
subplot(3,4,cont)
end

if (cont==5)   
subplot(3,4,cont-1)
end

if (cont<=9)&&(cont>5)   
subplot(3,4,cont-1)
end

if (cont==10)  
subplot(3,4,(cont-1)-1)
end

if (cont<=14)&&(cont>10)   
subplot(3,4,cont-2)
end

if (cont==15)   
subplot(3,4,cont-3)
end


imagesc(lon,lat,Y_Tropomi'); colormap jet;set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
clear lat lon ilat ilon tropomi j;
hold on

end


if (cont==4||cont==9||cont==14)
       sprintf('esta %i',i)
       
       if(cont==4)
       subplot(3,4,cont-1)
       end
       if(cont==9)
       subplot(3,4,cont-2)
       end
       if(cont==14)
       subplot(3,4,cont-3)
       end
       
       
       
       Y_Tropomi_2=Y_Tropomi;
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));
    
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{k,i-1,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
        end
   end 
   end
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_Tropomi(isnan(Y_Tropomi))=0;Y_Tropomi_2(isnan(Y_Tropomi_2))=0;
Sum_trop=Y_Tropomi+Y_Tropomi_2;Sum_trop(Sum_trop==0)=NaN;
 imagesc(lon,lat,(Sum_trop)'); colormap jet; % axis xy; 
set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
clear lat lon ilat ilon tropomi Y_Tropomi j;

end
   cont=cont+1; 

end



end

%%   Sigmas

figure
cont=1;
for k=1:3
for i=1:5
    
if ((cont~=4)&&(cont~=9)&&(cont~=14))
%     figure(i)
i
k
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'sigma'); %ask for sigma 
    
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));
    
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{k,i,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
        end
   end 
    end
if (cont<=4)   
subplot(3,4,cont)
end

if (cont==5)   
subplot(3,4,cont-1)
end

if (cont<=9)&&(cont>5)   
subplot(3,4,cont-1)
end

if (cont==10)  
subplot(3,4,(cont-1)-1)
end

if (cont<=14)&&(cont>10)   
subplot(3,4,cont-2)
end

if (cont==15)   
subplot(3,4,cont-3)
end


imagesc(lon,lat,Y_Tropomi'); colormap jet;set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,' Error instrument 1e15 mlc/cm2');caxis([0 3]);;title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
clear lat lon ilat ilon tropomi j;
hold on

end


if (cont==4||cont==9||cont==14)
       sprintf('esta %i',i)
       
       if(cont==4)
       subplot(3,4,cont-1)
       end
       if(cont==9)
       subplot(3,4,cont-2)
       end
       if(cont==14)
       subplot(3,4,cont-3)
       end
       
       
       
       Y_Tropomi_2=Y_Tropomi;
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));
    
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{k,i-1,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
        end
   end 
   end
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_Tropomi(isnan(Y_Tropomi))=0;Y_Tropomi_2(isnan(Y_Tropomi_2))=0;
Sum_trop=Y_Tropomi+Y_Tropomi_2;Sum_trop(Sum_trop==0)=NaN;
 imagesc(lon,lat,(Sum_trop)'); colormap jet; % axis xy; 
set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,' Error instrument 1e15 mlc/cm2');caxis([0 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
clear lat lon ilat ilon tropomi Y_Tropomi j;

end
   cont=cont+1; 

end



end