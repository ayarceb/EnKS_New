% code to plot the satellite files from the TROPOMI simulation of the
% EnKS-MC paper saving values of the random values extracted

% H_L   Tropomi without random extracted for assimilate

% H_V   Tropomi extracted to validate

% R_L   Sigma value without random extracted for assimilate

% R_V  Sigma extracted to validate

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


cont=1;
for k=1:3
for i=1:5
    
if ((cont~=4)&&(cont~=9)&&(cont~=14))
%     figure(i)
i
k
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
    sigma=ncread(files{i},'sigma');   % error satellite
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));  % H_L
    Y_Tropomi_sampled_validation=NaN(length(lon),length(lat));                                    % H_V
    Y_sigma=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));    %R_L
    Y_sigma_sampled_validation=NaN(length(lon),length(lat));                                       %R_V
       
    randam_obs{k,i,:}=randi([1 n],1,random_n(k));
   
    ccont=1;
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{k,i,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
      YY_H_L(cont,ccont)=tropomi(j);
      
      
      Y_sigma(ilon(j),ilat(j))=sigma(j); 
      YY_R_L(cont,ccont)=sigma(j);
      
      ccont=ccont+1;
       end
       
       
       if (ismember(j,randam_obs{k,i,:})==1)
      Y_Tropomi_sampled_validation(ilon(j),ilat(j))=tropomi(j); 
      YY_H_V=tropomi(j);
      
      Y_sigma_sampled_validation(ilon(j),ilat(j))=sigma(j); 
      YY_R_V=sigma(j);
      
      
       end
   end 
    end
    
    figure(5)
    imagesc(diag(YY_R_L(1,:)));colormap(flipud(hot));
    ylabel('Analysis states','FontSize',14);xlabel('Analysis states','FontSize',14);h=colorbar;ylabel(h,'Error instrument 1e15 mlc/cm2','FontSize',14)
    
    
%     figure(6)
%     for l=1:length(YY_R_L)
%         mean_R_L(1,l)=mean(YY_R_L);
%     end
%     
%     imagesc(diag(mean_R_L(1,:)));colormap(flipud(hot));
%     ylabel('Analysis states','FontSize',14);xlabel('Analysis states','FontSize',14);h=colorbar;ylabel(h,'Error instrument 1e15 mlc/cm2','FontSize',14)
%     caxis([0 1])
    
    figure(1)  
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
  
xlim([-76 -71.9]);ylim([9 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2','FontSize',10);caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
hold on

%%----

figure(2)

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


imagesc(lon,lat,Y_sigma'); colormap(flipud(hot));set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
xlim([-76 -71.9]);ylim([9 13]);h=colorbar;ylabel(h,'Error instruent 1e15 mlc/cm2','FontSize',10);caxis([-1 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')

hold on

%%-----

figure(3)

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


imagesc(lon,lat,Y_Tropomi_sampled_validation'); colormap jet;set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
xlim([-76 -71.9]);ylim([9 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2','FontSize',10);caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
hold on

%%-----

figure(4)

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


imagesc(lon,lat,Y_sigma_sampled_validation'); colormap(flipud(hot));set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
xlim([-76 -71.9]);ylim([9 13]);h=colorbar;ylabel(h,'Error instrument 1e15 mlc/cm2','FontSize',10);caxis([-1 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
clear lat lon ilat ilon tropomi j;
hold on


end

if (cont==4||cont==9||cont==14)
    sprintf('esta %i',i)
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
    sigma=ncread(files{i},'sigma');
          
   
    Y_Tropomi_2=Y_Tropomi;
    Y_Tropomi_sampled_validation_2=Y_Tropomi_sampled_validation;
    Y_sigma_2=Y_sigma;
    Y_sigma_sampled_validation_2=Y_sigma_sampled_validation;
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));  % H_L
    Y_Tropomi_sampled_validation=NaN(length(lon),length(lat));                                    % H_V
    Y_sigma=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));    %R_L
    Y_sigma_sampled_validation=NaN(length(lon),length(lat));                                       %R_V
       
    randam_obs{k,i,:}=randi([1 n],1,random_n(k));
  
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{k,i,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
      Y_sigma(ilon(j),ilat(j))=sigma(j); 
       end
       if (ismember(j,randam_obs{k,i,:})==1)
      Y_Tropomi_sampled_validation(ilon(j),ilat(j))=tropomi(j); 
      Y_sigma_sampled_validation(ilon(j),ilat(j))=sigma(j); 
      
       end
    end 
   end
   
    
     figure(1)
       if(cont==4)
       subplot(3,4,cont-1)
       end
       if(cont==9)
       subplot(3,4,cont-2)
       end
       if(cont==14)
       subplot(3,4,cont-3)
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
  
xlim([-76 -71.9]);ylim([9 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
   
    
     figure(2)
       if(cont==4)
       subplot(3,4,cont-1)
       end
       if(cont==9)
       subplot(3,4,cont-2)
       end
       if(cont==14)
       subplot(3,4,cont-3)
       end
    
    
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_sigma(isnan(Y_sigma))=0;Y_sigma_2(isnan(Y_sigma_2))=0;
Sum_trop_sigma=Y_sigma+Y_sigma_2;Sum_trop_sigma(Sum_trop_sigma==0)=NaN;
 imagesc(lon,lat,(Sum_trop_sigma)'); colormap(flipud(hot)); % axis xy; 
set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)

xlim([-76 -71.9]);ylim([9 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')


    
     figure(3)
       if(cont==4)
       subplot(3,4,cont-1)
       end
       if(cont==9)
       subplot(3,4,cont-2)
       end
       if(cont==14)
       subplot(3,4,cont-3)
       end
    
    
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_Tropomi_sampled_validation(isnan(Y_Tropomi_sampled_validation))=0;Y_Tropomi_sampled_validation_2(isnan(Y_Tropomi_sampled_validation_2))=0;
Sum_trop=Y_Tropomi_sampled_validation+Y_Tropomi_sampled_validation_2;Sum_trop(Sum_trop==0)=NaN;
imagesc(lon,lat,(Sum_trop)'); colormap jet; % axis xy; 
set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
xlim([-76 -71.9]);ylim([9 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
   
    
     figure(4)
       if(cont==4)
       subplot(3,4,cont-1)
       end
       if(cont==9)
       subplot(3,4,cont-2)
       end
       if(cont==14)
       subplot(3,4,cont-3)
       end
    
    
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_sigma_sampled_validation(isnan(Y_sigma_sampled_validation))=0;Y_sigma_sampled_validation_2(isnan(Y_sigma_sampled_validation_2))=0;
Sum_trop_sigma=Y_sigma_sampled_validation+Y_sigma_sampled_validation_2;Sum_trop_sigma(Sum_trop_sigma==0)=NaN;
 imagesc(lon,lat,(Sum_trop_sigma)'); colormap(flipud(hot)); % axis xy; 
set(gca,'YDir','normal'); hold on

S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)

xlim([-76 -71.9]);ylim([9 13]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')

clear lat lon ilat ilon tropomi Y_Tropomi j;hold on

end
   cont=cont+1; 

end

end
