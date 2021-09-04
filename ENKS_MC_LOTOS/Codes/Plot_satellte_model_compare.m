% code to plot the satellite files from the TROPOMI simulation of the
% EnKS-MC paper. Andrés Yarce Botero

clc;close all;clear all

addpath('/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC/ENKS_MC_LOTOS/Codes')
% Coordinates for the 59x63 array.
% Coordenadas        (latitude,longitude)
% Barranquilla       Latitud: 10.9878,    Longitud: -74.7889      (33,20)
% Santa Marta        Latitud: 11.24079,   Longitud: -74.19904      (37,26)
% Cartagena          Latitud: 10.3884731, Longitud: -75.488593     (27,12)
% Mina               Latitud: 9.5829857,  Longitud: -73.479376 ongitud: -74.7889      (33,20)
% Santa Marta        Latitud: 11.24079,   Longitud: -74.19904      (37,26)
% Cartagena          Latitud: 10.3884731, Longitud: -75.488593     (27,12)
% Mina               Latitud: 9.5829857,  Longitud: -73.479376      (18,35)
% Valledupar         Latitud: 10.4646054, Longitud: -73.2599952     ( 27,36)

ciudades={'Barranquilla ','Santa Marta ','Cartagena ','Mina Drummond ','Valledupar '};

lon_cities=[20,27,12,34,36];lat_cities=[33,37,27,18,27];

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

files_xb={'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190201_1800_xb.nc','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190202_1800_xb.nc',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1800_xb.nc','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1900_xb.nc',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190204_1900_xb.nc'};

% filesbase={'LE_Prueba_numero_4_EnKS_MC_column_20190201.nc','LE_Prueba_numero_4_EnKS_MC_column_20190202.nc',...
%     'LE_Prueba_numero_4_EnKS_MC_column_20190203.nc','LE_Prueba_numero_4_EnKS_MC_column_20190203.nc',...
%     'LE_Prueba_numero_4_EnKS_MC_column_20190204.nc'};

random_n=[1000,1500,2000];

cont=1;
k=1;
figure(1)
for i=1:5
    
if ((cont~=4))
     

    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');
    ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
    sigma=ncread(files{i},'sigma');
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);randam_obs=randi([1 n],1,random_n(k));
    Y_Tropomi_sampled_validation=NaN(length(lon),length(lat));
    
    Y_sigma=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));
    Y_sigma_sampled_validation=NaN(length(lon),length(lat));
    
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs)==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
       end
       if (ismember(j,randam_obs)==1)
      Y_Tropomi_sampled_validation(ilon(j),ilat(j))=tropomi(j); 
       end
    end
    
    end
    
    for j=1:length(ilat)
   
    if not(isnan(sigma(j)))
       if (ismember(j,randam_obs)==0)
      Y_sigma(ilon(j),ilat(j))=sigma(j); 
       end
       if (ismember(j,randam_obs)==1)
      Y_sigma_sampled_validation(ilon(j),ilat(j))=sigma(j); 
       end
   end 

    
    end
    
        
    Day_(:,:,i)=Y_Tropomi;
    Day_sampled_validation(:,:,i)=Y_Tropomi_sampled_validation;
    
    
    Day_sigma(:,:,i)=Y_sigma;
    Day_sampled_validation_sigma(:,:,i)=Y_sigma_sampled_validation;
    
if (cont<4)   
subplot(2,4,cont)
end

if (cont==5)   
subplot(2,4,cont-1)
end

lon_cities=[20,27,12,35,36];lat_cities=[33,37,27,18,27];

imagescwithnan(lon,lat,Y_Tropomi',jet,[0 0 0]); colormap jet;set(gca,'YDir','normal'); hold on
s=scatter(lon(lon_cities),lat(lat_cities),100,'s','MarkerEdgeColor',[1 1 1],...
              'MarkerFaceColor',[ 0.9961 0 0.9961],...
              'LineWidth',1.5,'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8); 

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1);
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0);
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °');
% saveas(fig,'Tropomi Retrieval','jpg')
% %..........................................................................................................................
%Second line subplots
if (cont<=4)   
subplot(2,4,cont+4)

end

if (cont==5)   
subplot(2,4,(cont-1)+4)
end


imagescwithnan(lon,lat,Y_Tropomi_sampled_validation',jet,[0 0 0]); ;xlim([-76 -71.5]);ylim([8 13]);set(gca,'YDir','normal'); hold on
s=scatter(lon(lon_cities),lat(lat_cities),100,'s','MarkerEdgeColor',[1 1 1],...
              'MarkerFaceColor',[ 0.9961 0 0.9961],...
              'LineWidth',1.5,'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8) 

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]),title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
xlim([-76 -71.5]);ylim([8 13])
clear ilat ilon tropomi j;
hold on

 end
% %.................................................................................................................
% %Day 3 two overpasses
% 
if (cont==4)
       sprintf('esta %i',i)
       
       if(cont==4)
       subplot(2,4,cont-1)
       end
              
       Y_Tropomi_2=Y_Tropomi;
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);randam_obs=randi([1 n],1,random_n(k));
    
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs)==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
        end
   end 
    end
%      
%     
%     
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_Tropomi(isnan(Y_Tropomi))=0;Y_Tropomi_2(isnan(Y_Tropomi_2))=0;
Sum_trop=Y_Tropomi+Y_Tropomi_2;Sum_trop(Sum_trop==0)=NaN;

    Day_(:,:,i)=Sum_trop;

imagesc(lon,lat,(Sum_trop)');  % axis xy; 
set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
s=scatter(lon(lon_cities),lat(lat_cities),100,'s','MarkerEdgeColor',[1 1 1],...
              'MarkerFaceColor',[ 0.9961 0 0.9961],...
              'LineWidth',1.5,'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8) 
          
% %-----------------------------------------------------------------------------------------------------------------------
%         
%    
if (cont==4)
       sprintf('esta %i',i)
       
       if(cont==4)
       subplot(2,4,cont-1+4)
       end
              
       Y_Tropomi_2_sampled_validation=Y_Tropomi_sampled_validation;
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
    
    Y_Tropomi_sampled_validation=NaN(length(lon),length(lat));n=length(ilat);randam_obs=randi([1 n],1,random_n(k));
    
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs)==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
       end
        if (ismember(j,randam_obs)==1)
      Y_Tropomi_sampled_validation(ilon(j),ilat(j))=tropomi(j); 
       end
   end 
    end
    
  
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_Tropomi_sampled_validation(isnan(Y_Tropomi_sampled_validation))=0;
Y_Tropomi_2_sampled_validation(isnan(Y_Tropomi_2_sampled_validation))=0;
Sum_trop_sampled_validation=Y_Tropomi_sampled_validation+Y_Tropomi_2_sampled_validation;
Sum_trop_sampled_validation(Sum_trop_sampled_validation==0)=NaN;
Day_sampled_validation(:,:,i)=Sum_trop_sampled_validation; 

imagescwithnan(lon,lat,Sum_trop_sampled_validation',jet,[0 0 0]);  % axis xy; 
set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)

xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
s=scatter(lon(lon_cities),lat(lat_cities),100,'s','MarkerEdgeColor',[1 1 1],...
              'MarkerFaceColor',[ 0.9961 0 0.9961],...
              'LineWidth',1.5,'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8) 
          
% %-----------------------------------------------------------------------------------------------------------------------
% 
 end       
%           
% clear ilat ilon tropomi Y_Tropomi j;
% 
 end
% 
 cont=cont+1; 
% 
 end
% 
% 
% Esta parte siguiente genera la visualización similar para los dias del 1 al
% 4 con sus valores del muestreo
% %% Plot values of the observation through time
% 
figure(2)
% 
for i=1:5
if(i<4)    
subplot(2,4,i)
imagescwithnan(lon,lat,squeeze(Day_(:,:,i)'),jet,[0 0 0]);axis xy;
xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
subplot(2,4,i+4)
imagescwithnan(lon,lat,squeeze(Day_sampled_validation(:,:,i)'),jet,[0 0 0]);axis xy;
xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')



Data_Sat_Assim(:,:,i)=Day_(:,:,i)
Data_Sat_Assim_validation(:,:,i)=Day_sampled_validation(:,:,i)

end


if(i==4)    
subplot(2,4,i-1)
imagescwithnan(lon,lat,squeeze(Day_(:,:,i)'),jet,[0 0 0]);axis xy;
xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
subplot(2,4,i+3)
imagescwithnan(lon,lat,squeeze(Day_sampled_validation(:,:,i)'),jet,[0 0 0]);axis xy;
xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')



Data_Sat_Assim(:,:,i)=Day_(:,:,i)
Data_Sat_Assim_validation(:,:,i)=Day_sampled_validation(:,:,i)

end



if(i>4)
subplot(2,4,i-1)
imagescwithnan(lon,lat,squeeze(Day_(:,:,i)'),jet,[0 0 0]);axis xy;
xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
subplot(2,4,i+3)
imagescwithnan(lon,lat,squeeze(Day_sampled_validation(:,:,i)'),jet,[0 0 0]);axis xy;
xlim([-76 -71.5]);ylim([8 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')

    
    Data_Sat_Assim(:,:,i-1)=Day_(:,:,i)
Data_Sat_Assim_validation(:,:,i-1)=Day_sampled_validation(:,:,i)

 end

 end
% 
% 


%%
figure

sz=20;
for j=1:5
    
    if(j<4)
        
  subplot(2,5,j)
    y=ncread(files_xb{j},'y'); % modelo
    yr=ncread(files{j},'yr');  % tropomi
%  scatter(y,yr,sz,'MarkerEdgeColor',[.5 .5 .5],...
%               'MarkerFaceColor',[0.1 0.1 0.1],...
%               'LineWidth',.5);
y(isnan(y)) = [];
yr(isnan(yr)) = [];
tbl = table(y',yr');
mdl = fitlm(tbl,'Var1 ~ Var2')
if(j~=3)
plotAdded(mdl)
 
hold on
  xlabel('Simulated retrieval mlc/cm2');
  ylabel('Tropospheric Column Density TROPOMI mlc/cm2');
  xlim([0 8])
  ylim([0 8])
  x = [0 8];
yy = [0 8];
% line(x,yy,'Color','g','LineStyle','--')
title(date{j});
subplot(2,5,5+j)
plot(y-yr,'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
xlabel('pixels')
  x = [0 3000];
yy = [0 0];
line(x,yy,'Color','g','LineStyle','--')

hold on
ylabel('mlc/cm2')
title('innovation (Y-H(x))')
grid on
 end      
    end
    
    
     if(j==4)
         y_3=y;
         yr_3=yr;
  subplot(2,5,j-1)
    y=[ncread(files_xb{j},'y') y_3]; % modelo
    yr=[ncread(files{j},'yr') yr_3];  % tropomi
%      scatter(y,yr,sz,'MarkerEdgeColor',[.5 .5 .5],...
%               'MarkerFaceColor',[0.1 0.1 0.1],...
%               'LineWidth',.5);
y(isnan(y)) = [];
yr(isnan(yr)) = [];
tbl = table(y',yr');
mdl = fitlm(tbl,'Var1 ~ Var2')
plotAdded(mdl)

  xlabel('Simulated retrieval mlc/cm2');
  ylabel('Tropospheric Column Density TROPOMI mlc/cm2');
  xlim([0 8])
  ylim([0 8])
  x = [0 8];
yy = [0 8];
% line(x,yy,'Color','g','LineStyle','--')
title(date{j});
subplot(2,5,4+j)
plot(y-yr,'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
xlabel('pixels')
  x = [0 3000];
yy = [0 0];
line(x,yy,'Color','g','LineStyle','--')

ylabel('mlc/cm2')
title('innovation (Y-H(x))')
grid on
     end
    
     if(j>4)
  subplot(2,5,j-1)
    y=ncread(files_xb{j},'y'); % modelo
    yr=ncread(files{j},'yr');  % tropomi
%      scatter(y,yr,sz,'MarkerEdgeColor',[.5 .5 .5],...
%               'MarkerFaceColor',[0.1 0.1 0.1],...
%               'LineWidth',.5);
y(isnan(y)) = [];
yr(isnan(yr)) = [];
tbl = table(y',yr');
mdl = fitlm(tbl,'Var1 ~ Var2')
plotAdded(mdl,'ro')
  xlabel('Simulated retrieval mlc/cm2');
  ylabel('Tropospheric Column Density TROPOMI mlc/cm2');
  xlim([0 8])
  ylim([0 8])
  x = [0 8];
yy = [0 8];
% line(x,yy,'Color','g','LineStyle','--')
title(date{j});
subplot(2,5,4+j)
plot(y-yr,'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
xlabel('pixels')
  x = [0 3000];
yy = [0 0];
% line(x,yy,'Color','g','LineStyle','--')

ylabel('mlc/cm2')
title('innovation (Y-H(x))')
grid on
    end
     
     
    end






 


