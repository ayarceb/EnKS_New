%% Code to plot the satellite files from the TROPOMI simulation of the
% EnKS-MC paper saving values of the random values extracted
% Andres Yarce Botero  Santiago Lopez 2021

% This code saves the following matrixes for have for the analysis
% assimilation procedure:

% H_L_TROPOMI(i,:,:) Tropomi without random extracted for assimilate for each day
% R_L_TROPOMI(i,:,:) Sigma value without random extracted for assimilate for each day
% H_V_TROPOMI(i,:,:) Tropomi extracted to validate for each day
% R_V_TROPOMI(i,:,:) Sigma extracted to validate for each day

%------------------------------------------------------------------------------------------------------------------------------
% Variable characteristics    Size no2:       58x63x1x97
% From the LOTOS-EUROS        Dimensions: longitude,latitude,level,time
%                             Datatype:   single
%                             Attributes:
%                             Standard_name   = 'mole_fraction_of_nitrogen_dioxide_in_air'
%                             Long_name       = 'volume mixing ratio of NO2 in humid air'
%                             Units           = 'mole mole-1'
%                             CoordinateAxes = 'time level latitude longitude'
%                             Molemass        = 0.046005
%                             Molemass_unit   = 'kg mole-1'
% Coordinates for the 58x63 array.
% Coordenadas        (latitude,longitude)
% Barranquilla       Latitud: 10.9878,    Longitud: -74.7889     (33,20)
% Santa Marta        Latitud: 11.24079,   Longitud: -74.19904    (37,26)
% Cartagena          Latitud: 10.3884731, Longitud: -75.488593   (27,12)
% Mina               Latitud: 9.5829857,  Longitud: -73.479376   (18,35)
% Valledupar         Latitud: 10.4646054, Longitud: -73.2599952  ( 27,36)
%--------------------------------------------------------------------------------------------------------------------------------

clc;close all;clear all
addpath('/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC_new')
addpath('/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC_new/ENKS_MC_LOTOS/Codes')
date={'2019-02-01','2019-02-02','2019-02-03','2019-02-03','2019-02-04'};
tit={'NO_2 LOTOS-EUROS Column','Vertical column density','Simulated retrieval'};
ciudades={'Barranquilla ','Santa Marta ','Cartagena ','Mina Drummond ','Valledupar '};
lon=[20,26,12,35,36];lat=[33,37,27,18,27];  % index in this matrix of the cities
mydir='/run/media/dirac/Datos/scratch/projects/Prueba_numero_4_EnKS_MC/Prueba_numero_4_EnKS_MC/output';
cd(mydir)

files={'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190201_1800.nc','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190202_1800.nc',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1800.nc','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1900.nc',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190204_1900.nc'};

filesbase={'LE_Prueba_numero_4_EnKS_MC_column_20190201_xb.nc','LE_Prueba_numero_4_EnKS_MC_column_20190202_xb.nc',...
    'LE_Prueba_numero_4_EnKS_MC_column_20190203_xb.nc','LE_Prueba_numero_4_EnKS_MC_column_20190203_xb.nc',...
    'LE_Prueba_numero_4_EnKS_MC_column_20190204_xb.nc'};

files_xb={'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190201_1800_xb.nc','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190202_1800_xb.nc',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1800_xb.nc','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1900_xb.nc',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190204_1900_xb.nc'};

%%---------------------------------------------------------------------------------------------------------------------------
% estos siguientes path están comentados debido a que es la corrida Prueba
% 6 mayor latitud:

% mydir='/run/media/dirac/Datos/scratch/projects/Prueba_numero_6_EnKS_MC/Prueba_numero_6_EnKS_MC/output';
% cd(mydir)
% 
% files={'LE_Prueba_numero_6_EnKS_MC_tropomi-no2_20190201_1800.nc','LE_Prueba_numero_6_EnKS_MC_tropomi-no2_20190202_1800.nc',...
%     'LE_Prueba_numero_6_EnKS_MC_tropomi-no2_20190203_1800.nc','LE_Prueba_numero_6_EnKS_MC_tropomi-no2_20190203_1900.nc',...
%     'LE_Prueba_numero_6_EnKS_MC_tropomi-no2_20190204_1900.nc'};
% 
% filesbase={'LE_Prueba_numero_6_EnKS_MC_column_20190201_xb.nc','LE_Prueba_numero_6_EnKS_MC_column_20190202_xb.nc',...
%     'LE_Prueba_numero_6_EnKS_MC_column_20190203_xb.nc','LE_Prueba_numero_6_EnKS_MC_column_20190203_xb.nc',...
%     'LE_Prueba_numero_6_EnKS_MC_column_20190204_xb.nc'};
%%---------------------------------------------------------------------------------------------------------------------------
%% Part for import the TROPOMI retrieval and the associated error sigma value to sampled based on a number of samples
k=1;
cont=1;
PR=0.9;   % Adjust PR from 0 to 1 percentage sampled data for validation 
for i=1:5  % son 4 dias pero el tercer dia tiene 2 overpasses
    
if (cont~=4)  %This if collect the days different than cont=4
    
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude');
    ilat=ncread(files{i},'ilat');ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');    % Tropospheric vertical column of nitrogen dioxide
    sigma=ncread(files{i},'sigma');   % error satellite
    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);                                        % H_L
    Y_Tropomi_sampled_validation=NaN(length(lon),length(lat));                                    % H_V
    Y_sigma=NaN(length(lon),length(lat));n=length(ilat);                                          %R_L
    Y_sigma_sampled_validation=NaN(length(lon),length(lat));                                      %R_V
       
    randam_obs{i,:}=randperm(n,floor(n*PR));  % almacenar las posiciones de las posiciones random para cada uno de los valores 
   
    ccont=1;
    
 if(i==1)||(i==2)||(i==5) 
      V=cell2mat(randam_obs(i)); 
      SS{i,:}=sigma(V);   % SS --> Sigma error satellite
      TT{i,:}=tropomi(V); % TT --> Satellite value
 for j=1:length(ilat)
  
     % the following two conditional if extract the value from the
     % assimilation and save this values for validation step
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{i,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
%       YY_H_L(cont,ccont)=tropomi(j);
      Y_sigma(ilon(j),ilat(j))=sigma(j); 
%       YY_R_L(cont,ccont)=sigma(j);
%       ccont=ccont+1;
       end
         
       if (ismember(j,randam_obs{i,:})==1)
      Y_Tropomi_sampled_validation(ilon(j),ilat(j))=tropomi(j); 
%       YY_H_V(cont,c_cont)=tropomi(j);
      Y_sigma_sampled_validation(ilon(j),ilat(j))=sigma(j); 
%       YY_R_V(cont,c_cont)=sigma(j);
%       c_cont=c_cont+1;
   end
   end 
 end
 end
 
 if(i==3)    
     V=cell2mat(randam_obs(cont)); 
     SS{i,:}=sigma(V);
     TT{i,:}=tropomi(V); % TT --> Satellite value
 for j=1:length(ilat)
  
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{i,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
      Y_sigma(ilon(j),ilat(j))=sigma(j); 

      
      ccont=ccont+1;
       end
              
       if (ismember(j,randam_obs{i,:})==1)
      Y_Tropomi_sampled_validation(ilon(j),ilat(j))=tropomi(j); 
      Y_sigma_sampled_validation(ilon(j),ilat(j))=sigma(j); 

       
   end
   end 
 end
 end
  
%% Section to plot the Covariance observartion matrix 
%  figure(5)
%     imagesc(diag(YY_R_L(1,:)));colormap(flipud(hot));
%     ylabel('Analysis states','FontSize',14);xlabel('Analysis states','FontSize',14);h=colorbar;ylabel(h,'Error instrument 1e15 mlc/cm2','FontSize',14)
%     
%    
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
subplot(1,4,cont)
end

if (cont==5)   
subplot(1,4,cont-1)
end

H_L_TROPOMI(i,:,:)=Y_Tropomi';

imagesc(lon,lat,Y_Tropomi'); colormap jet;set(gca,'YDir','normal'); hold on

S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
xlim([-76.5 -71.28]);ylim([9 13.65]);h=colorbar;ylabel(h,'1e15 mlc/cm2','FontSize',10);caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
hold on
%%----
figure(2)

if (cont<=4)   
subplot(1,4,cont)
end

if (cont==5)   
subplot(1,4,cont-1)
end

R_L_TROPOMI(i,:,:)=Y_sigma';
imagesc(lon,lat,Y_sigma'); colormap(flipud(hot));set(gca,'YDir','normal'); hold on

S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
xlim([-76.5 -71.28]);ylim([9 13.65]);h=colorbar;ylabel(h,'Error instruent 1e15 mlc/cm2','FontSize',10);caxis([-1 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')

hold on
%%-----
figure(3)

if (cont<=4)   
subplot(1,4,cont)
end

if (cont==5)   
subplot(1,4,cont-1)
end

H_V_TROPOMI(i,:,:)=Y_Tropomi_sampled_validation';
imagesc(lon,lat,Y_Tropomi_sampled_validation'); colormap jet;set(gca,'YDir','normal'); hold on

S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
xlim([-76.5 -71.28]);ylim([9 13.65]);h=colorbar;ylabel(h,'1e15 mlc/cm2','FontSize',10);caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
hold on
%%-----
figure(4)

if (cont<=4)   
subplot(1,4,cont)
end

if (cont==5)   
subplot(1,4,cont-1)
end

R_V_TROPOMI(i,:,:)=Y_sigma_sampled_validation';
imagesc(lon,lat,Y_sigma_sampled_validation'); colormap(flipud(hot));set(gca,'YDir','normal'); hold on
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
xlim([-76.5 -71.28]);ylim([9 13.65]);h=colorbar;ylabel(h,'Error instrument 1e15 mlc/cm2','FontSize',10);caxis([-1 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
clear lat lon ilat ilon tropomi j;
hold on

end

% Esta segunda parte es el "Pain in the ass" ya que el tercer dia de datos
% tiene dos overpasses, por lo que hay que a tomar el del dia pasado y
% adicionarlo al pedazo actual

if (cont==4)
   
    sprintf('esta %i',i)
    lat=ncread(filesbase{i},'latitude');lon=ncread(filesbase{i},'longitude'); ilat=ncread(files{i},'ilat');
    ilon=ncread(files{i},'ilon');
    tropomi=ncread(files{i},'yr');
    sigma=ncread(files{i},'sigma');n=length(ilat)
    randam_obs{i,:}=randperm(n,floor(n*PR));  % almacenar las posiciones de las posiciones random para cada uno de los valores 
    V=cell2mat(randam_obs(i)); 
    SS{i,:}=sigma(V);   % SS --> Sigma error satellite
    TT{i,:}=tropomi(V); % TT --> Satellite value
%---se almacenan los pasados inicialmente--------------------------
    Y_Tropomi_2=Y_Tropomi; 
    Y_Tropomi_sampled_validation_2=Y_Tropomi_sampled_validation;
    Y_sigma_2=Y_sigma;
    Y_sigma_sampled_validation_2=Y_sigma_sampled_validation;
%------------------------------------------------------------------    
    Y_Tropomi=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));  % H_L
    Y_Tropomi_sampled_validation=NaN(length(lon),length(lat));                                    % H_V
    Y_sigma=NaN(length(lon),length(lat));n=length(ilat);%randam_obs=randi([1 n],1,random_n(k));    %R_L
    Y_sigma_sampled_validation=NaN(length(lon),length(lat));                                       %R_V
       
 
    for j=1:length(ilat)
   
    if not(isnan(tropomi(j)))
       if (ismember(j,randam_obs{i-1,:})==0)
      Y_Tropomi(ilon(j),ilat(j))=tropomi(j); 
      Y_sigma(ilon(j),ilat(j))=sigma(j); 
       end
       if (ismember(j,randam_obs{i-1,:})==1)
      Y_Tropomi_sampled_validation(ilon(j),ilat(j))=tropomi(j); 
      Y_sigma_sampled_validation(ilon(j),ilat(j))=sigma(j); 
      
       end
    end 
   end
    
     figure(1)
       if(cont==4)
       subplot(1,4,cont-1)
       end
    
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_Tropomi(isnan(Y_Tropomi))=0;Y_Tropomi_2(isnan(Y_Tropomi_2))=0;
Sum_trop=Y_Tropomi+Y_Tropomi_2;Sum_trop(Sum_trop==0)=NaN;

H_L_TROPOMI(i,:,:)=Sum_trop';
imagesc(lon,lat,(Sum_trop)'); colormap jet; % axis xy; 
set(gca,'YDir','normal'); hold on

S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
xlim([-76.5 -71.28]);ylim([9 13.65]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
   
    
     figure(2)
       if(cont==4)
       subplot(1,4,cont-1)
       end
         
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_sigma(isnan(Y_sigma))=0;Y_sigma_2(isnan(Y_sigma_2))=0;
Sum_trop_sigma=Y_sigma+Y_sigma_2;Sum_trop_sigma(Sum_trop_sigma==0)=NaN;
R_L_TROPOMI(i,:,:)=Sum_trop_sigma';
imagesc(lon,lat,(Sum_trop_sigma)'); colormap(flipud(hot)); % axis xy; 
set(gca,'YDir','normal'); hold on

S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)

xlim([-76.5 -71.28]);ylim([9 13.65]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')
    
     figure(3)
       if(cont==4)
       subplot(1,4,cont-1)
       end
        
    
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_Tropomi_sampled_validation(isnan(Y_Tropomi_sampled_validation))=0;Y_Tropomi_sampled_validation_2(isnan(Y_Tropomi_sampled_validation_2))=0;
Sum_trop=Y_Tropomi_sampled_validation+Y_Tropomi_sampled_validation_2;Sum_trop(Sum_trop==0)=NaN;


H_V_TROPOMI(i,:,:)=Sum_trop';
imagesc(lon,lat,(Sum_trop)'); colormap jet; % axis xy; 
set(gca,'YDir','normal'); hold on

S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
  
xlim([-76.5 -71.28]);ylim([9 13.65]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 7]);title(date{i});xlabel('longitude °');ylabel('latitude °')
   
    
     figure(4)
       if(cont==4)
       subplot(1,4,cont-1)
       end
       
    
%  imagesc(lon,lat,(Y_Tropomi_2)');hold on;
Y_sigma_sampled_validation(isnan(Y_sigma_sampled_validation))=0;Y_sigma_sampled_validation_2(isnan(Y_sigma_sampled_validation_2))=0;
Sum_trop_sigma=Y_sigma_sampled_validation+Y_sigma_sampled_validation_2;Sum_trop_sigma(Sum_trop_sigma==0)=NaN;
 
R_V_TROPOMI(i,:,:)=Sum_trop_sigma';
imagesc(lon,lat,(Sum_trop_sigma)'); colormap(flipud(hot)); % axis xy; 
set(gca,'YDir','normal'); hold on

S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)

xlim([-76.5 -71.28]);ylim([9 13.65]);h=colorbar;ylabel(h,'1e15 mlc/cm2');caxis([-1 3]);title(date{i});xlabel('longitude °');ylabel('latitude °')

clear lat lon ilat ilon tropomi Y_Tropomi j;hold on

end
   cont=cont+1; 

end

save('/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/Randam.mat','randam_obs')
save('/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/Sigmas.mat','SS')
save('/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/Tropomi.mat','TT')
save('/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/Matrices.mat','H_L_TROPOMI','R_L_TROPOMI','H_V_TROPOMI','R_V_TROPOMI')
%% Plot part to corroborate that the Observation and sampling observation operator write appropiately
% For the next matrices the time 4 are the sum of two overpasses
% close all
% for i=1:5
% figure(7)
% subplot(1,5,i); imagesc(squeeze(H_L_TROPOMI(i,:,:)));colorbar; axis xy
% sgtitle('H L (For assimilation)')
% figure(8)
% subplot(1,5,i); imagesc(squeeze(R_L_TROPOMI(i,:,:)));colorbar;axis xy
% sgtitle('R L (For assimilation)')
% figure(9)
% subplot(1,5,i); imagesc(squeeze(H_V_TROPOMI(i,:,:)));colorbar; axis xy
% sgtitle('H V (For validation)')
% figure(10)
% subplot(1,5,i); imagesc(squeeze(R_V_TROPOMI(i,:,:)));colorbar;axis xy
% sgtitle('H V (For validation)')
% end
%% Import model and tropomi from satellite files to calculate Innovations (Y-H(x)). Import files from the LOTOS-EUROS simulation and 
 close all;
figure(10)

for j=1:5
    V=cell2mat(randam_obs(j));
    
    if((j==1)||(j==2))
        
    subplot(3,5,j)
    y=ncread(files_xb{j},'y'); % modelo
    yr=ncread(files{j},'yr');  % tropomi
    y_model=y;
    yr_tropomi=yr;
    y(randam_obs{j,:})=0;
    yr(randam_obs{j,:})=0;
    %  scatter(y,yr,sz,'MarkerEdgeColor',[.5 .5 .5],...
    %               'MarkerFaceColor',[0.1 0.1 0.1],...
    %               'LineWidth',.5);
    % y(isnan(y)) = [];
    % yr(isnan(yr)) = [];
    tbl = table(y',yr');
    mdl = fitlm(tbl,'Var1 ~ Var2')
    
    
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
    subplot(3,5,5+j)
    d{j,:}=y-yr;  % Innovation
    plot(squeeze(d{j,:}),'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
           
    d{j,:}(d{j,:}==0)=[];
    xlabel('pixels')
    x = [0 3000];
    yy = [0 0];
    % line(x,yy,'Color','g','LineStyle','--')
    hold on;ylabel('mlc/cm2');title('innovation (Y-H(x))');grid on


    subplot(3,5,10+j)   
    d_val{j,:}=y_model-yr_tropomi;
    A=d_val{j,:}(V);
    d_valid{j,1}=A;  % Innovation_validation
    plot(d_val{j,:},'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
    xlabel('pixels')
    x = [0 3000];
    yy = [0 0];
    % line(x,yy,'Color','g','LineStyle','--')
    hold on;ylabel('mlc/cm2');title('innovation (Y-H(x)) without sampled');grid on

    end   


    if(j==3)
        
    subplot(3,5,j)
    y=ncread(files_xb{j},'y'); % modelo
    yr=ncread(files{j},'yr');  % tropomi
    y_model=y;
    yr_tropomi=yr;
    y(randam_obs{j,:})=0;
    yr(randam_obs{j,:})=0;
    %  scatter(y,yr,sz,'MarkerEdgeColor',[.5 .5 .5],...
    %               'MarkerFaceColor',[0.1 0.1 0.1],...
    %               'LineWidth',.5);
    % y(isnan(y)) = [];
    % yr(isnan(yr)) = [];
    tbl = table(y',yr');
    mdl = fitlm(tbl,'Var1 ~ Var2')
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


    subplot(3,5,5+j)
    d{j,:}=y-yr;   %Innovation
    plot(squeeze(d{j,:}),'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
    d{j,:}(d{j,:}==0)=[];
    xlabel('pixels')
    x = [0 3000];
    yy = [0 0];
    % line(x,yy,'Color','g','LineStyle','--')
    hold on;ylabel('mlc/cm2');title('innovation (Y-H(x))');grid on


    subplot(3,5,10+j)
    
    
    d_val{j,:}=y_model-yr_tropomi;
    A=d_val{j,:}(V);
    d_valid{j,1}=A;  % Innovation_validation
    plot(d_val{j,:},'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
    xlabel('pixels')
    x = [0 3000];
    yy = [0 0];
    % line(x,yy,'Color','g','LineStyle','--')
    hold on;ylabel('mlc/cm2');title('innovation (Y-H(x)) without sampled');grid on

    end 

   
   
if(j==4)
%     y_3=y;
%     yr_3=yr;
    subplot(3,5,j)
    y=[ncread(files_xb{j},'y')]; % modelo
    yr=[ncread(files{j},'yr')];  % tropomi
    y_model=y;
    yr_tropomi=yr;
    y(randam_obs{j,:})=0;
    yr(randam_obs{j,:})=0;
    %      scatter(y,yr,sz,'MarkerEdgeColor',[.5 .5 .5],...
    %               'MarkerFaceColor',[0.1 0.1 0.1],...
    %               'LineWidth',.5);
%     y(isnan(y)) = [];
%     yr(isnan(yr)) = [];
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
    subplot(3,5,5+j)
    d{j,:}=y-yr;    %Innovation
    plot(squeeze(d{j,:}),'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
   d{j,:}(d{j,:}==0)=[];
   xlabel('pixels')
   x = [0 3000];yy = [0 0];line(x,yy,'Color','g','LineStyle','--')

   ylabel('mlc/cm2')
   title('innovation (Y-H(x))')
   grid on
   
   subplot(3,5,10+j)
    d_val{j,:}=y_model-yr_tropomi;
    A=d_val{j,:}(V);
    d_valid{j,1}=A;  % Innovation_validation
    plot(d_val{j,:},'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
    xlabel('pixels')
    x = [0 3000];
    yy = [0 0];
    % line(x,yy,'Color','g','LineStyle','--')
    hold on;ylabel('mlc/cm2');title('innovation (Y-H(x)) without sampled');grid on
end
    
    if(j>4)
    subplot(3,5,5)
    y=ncread(files_xb{j},'y'); % modelo
    yr=ncread(files{j},'yr');  % tropomi
    %      scatter(y,yr,sz,'MarkerEdgeColor',[.5 .5 .5],...
    %               'MarkerFaceColor',[0.1 0.1 0.1],...
    %               'LineWidth',.5);
    y_model=y;
    yr_tropomi=yr;
    y(randam_obs{j,:})=0;
    yr(randam_obs{j,:})=0;
%   y(isnan(y)) = [];
%   yr(isnan(yr)) = [];
    tbl = table(y',yr');
    mdl = fitlm(tbl,'Var1 ~ Var2')
    plotAdded(mdl)
    xlabel('Simulated retrieval mlc/cm2');
    ylabel('Tropospheric Column Density TROPOMI mlc/cm2');
    xlim([0 8]);ylim([0 8])
    x = [0 8];yy = [0 8];
    % line(x,yy,'Color','g','LineStyle','--')
    title(date{j});
    
    subplot(3,5,5+j)
    d{j,:}=y-yr;
    plot(squeeze(d{j,:}),'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
    d{j,:}(d{j,:}==0)=[];
    xlabel('pixels')
        x = [0 3000];
    yy = [0 0];
    ylabel('mlc/cm2');title('innovation (Y-H(x))');grid on
  
    subplot(3,5,10+j)
    d_val{j,:}=y_model-yr_tropomi;
    A=d_val{j,:}(V);
    d_valid{j,1}=A;  % Innovation_validation
    plot(d_val{j,:},'ro','MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[0.1 0.1 0.1],...
               'LineWidth',.001,'MarkerSize',2)
    xlabel('pixels')
    x = [0 3000];
    yy = [0 0];
    % line(x,yy,'Color','g','LineStyle','--')
    hold on;ylabel('mlc/cm2');title('innovation (Y-H(x)) without sampled');grid on
    end   
end

save('/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/innovation.mat','d','d_valid','d_val')

%###########################################################################################################################
%% From code EnKS Assimilate  part for time plots and cholesky decomposition as well reshapes to 2D 
%###########################################################################################################################
% 
ciudades={'Barranquilla ','Santa Marta ','Cartagena ','Mina Drummond ','Valledupar '};

lon=[20,26,12,35,36];lat=[33,37,27,18,27];

name_run='Prueba_numero_4_EnKS_MC';

mydir='/run/media/dirac/Datos/scratch/projects/Prueba_numero_4_EnKS_MC/Prueba_numero_4_EnKS_MC/output';
cd(mydir)

lati=ncread('LE_Prueba_numero_4_EnKS_MC_column_20190205_xi02a.nc','latitude');
long=ncread('LE_Prueba_numero_4_EnKS_MC_column_20190205_xi02a.nc','longitude');

mydir='/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC_new/ENKS_MC_LOTOS/Codes';
cd(mydir)

%%===Path where the output are located from the EnKS_MC ensemble first propagation
mydir=append('/run/media/dirac/Datos/scratch/projects/',name_run,'/',name_run,'/output');cd(mydir)

system('rm no2_column_ens*.nc');system('rm Merge_*.nc'); system('rm Ens_dc_*.nc');
cd ..
system('mv LE_Prueba_numero_4_EnKS_MC_dc_*.nc output')
mydir=append('/run/media/dirac/Datos/scratch/projects/',name_run,'/',name_run,'/output');cd(mydir)

no2=zeros(58,63,97,40); no2_dc=zeros(58,63,97,40);

t1 = datetime(2019,2,1,0,0,0);t2 = datetime(2019,2,5,0,0,0);t = t1:hours(1):t2;
figure
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

mydir=append('/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC_new/ENKS_MC_LOTOS');cd(mydir)

system('./Merging_dc.sh')  % Bash code to merge dc factors

mydir=append('/run/media/dirac/Datos/scratch/projects/',name_run,'/',name_run,'/output');cd(mydir)

%% Take state snapshots of the states from the forward model

NO2_state_0=squeeze(NO2_state(:,:,1,:));
NO2_state_1=squeeze(NO2_state(:,:,17,:));
NO2_state_2=squeeze(NO2_state(:,:,41,:));
NO2_state_3=squeeze(NO2_state(:,:,65,:));
NO2_state_4=squeeze(NO2_state(:,:,97,:));

NNO2_state_0=reshape(NO2_state_0,58*63,40);
NNO2_state_1=reshape(NO2_state_1,58*63,40);
NNO2_state_2=reshape(NO2_state_2,58*63,40);
NNO2_state_3=reshape(NO2_state_3,58*63,40);
NNO2_state_4=reshape(NO2_state_4,58*63,40);

%% Calculate Cholesky for the states, Plot Cholesky and Square Cholesky  (Only states!)
% figure
STATES_lag0=[NNO2_state_0];
STATES_lag1=[NNO2_state_1];
STATES_lag2=[NNO2_state_2];
STATES_lag3=[NNO2_state_3];
STATES_lag4=[NNO2_state_4];

% B_lag0=Calculo_B_Cholesky(STATES_lag0,1);
% B_lag1=Calculo_B_Cholesky(STATES_lag1,1);
% B_lag2=Calculo_B_Cholesky(STATES_lag2,1);
% B_lag3=Calculo_B_Cholesky(STATES_lag3,1);
%B_lag4=Calculo_B_Cholesky(STATES_lag4,1);

[B_chol0,B_square0]=Calculo_B_Cholesky(STATES_lag0,1);
[B_chol1,B_square1]=Calculo_B_Cholesky(STATES_lag1,1);
[B_chol2,B_square2]=Calculo_B_Cholesky(STATES_lag2,1);
[B_chol3,B_square3]=Calculo_B_Cholesky(STATES_lag3,1);
[B_chol4,B_square4]=Calculo_B_Cholesky(STATES_lag4,1);


save('/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/B_square.mat','B_square0','B_square1','B_square2','B_square3','B_square4')
figure(20)
subplot(1,3,2) 
imagesc(B_square)
colorbar;colormap(jet)
caxis([0 1e-9]);ylabel('States');xlabel('States')
title('Square of the Cholesky matrix')
subplot(1,3,3)
imagesc(B_chol)
colorbar;colormap(jet)
caxis([0 1e-16]);ylabel('States');xlabel('States')
title('Complete Cholesky matrix')


%%  Take state snapshots of the states from the forward model PLOT DC

figure

for j=1:5
for  n=1:ens 
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
% mydir=append('/run/media/dirac/Datos/scratch/projects/',name_run,'/',name_run,'/output');cd(mydir)

% %% Take snapshots of the DC factor from the forward model
% 
% % 
% % NNO2_dc_0=reshape(NO2_dc_0,58*63,40);
% % NNO2_dc_1=reshape(NO2_dc_1,58*63,40);
% % NNO2_dc_2=reshape(NO2_dc_2,58*63,40);
% % NNO2_dc_3=reshape(NO2_dc_3,58*63,40);
% % NNO2_dc_4=reshape(NO2_dc_4,58*63,40);
% % 
% % figure
% % dc_lag0=[NNO2_dc_0(1,:)];
% % dc_lag1=[NNO2_dc_1(1,:)];
% % dc_lag2=[NNO2_dc_2(1,:)];
% % dc_lag3=[NNO2_dc_3(1,:)];
% % dc_lag4=[NNO2_dc_4(1,:)];
% % 
% % B_dc_lag0=Calculo_B_Cholesky(dc_lag0,1);
% % B_dc_lag1=Calculo_B_Cholesky(dc_lag1,1);
% % B_dc_lag2=Calculo_B_Cholesky(dc_lag2,1);
% % B_dc_lag3=Calculo_B_Cholesky(dc_lag3,1);
% % B_dc_lag4=Calculo_B_Cholesky(dc_lag4,1);
% % 
% % subplot(2,3,1);imagesc(B_dc_lag0);colorbar;colormap(jet);title('L 0');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% % subplot(2,3,2);imagesc(B_dc_lag1);colorbar;colormap(jet);title('L 1');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% % subplot(2,3,3);imagesc(B_dc_lag2);colorbar;colormap(jet);title('L 2');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% % subplot(2,3,4);imagesc(B_dc_lag3);colorbar;colormap(jet);title('L 3');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% % subplot(2,3,5);imagesc(B_dc_lag4);colorbar;colormap(jet);title('L 4');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% % States + parameteres and calculate EnKF background covariance matrix ensemble approximation
% t_index=97;
% 
% NO2_dc_0=squeeze(no2_dc(:,:,1,:));
% NO2_dc_1=squeeze(no2_dc(:,:,17,:));
% NO2_dc_2=squeeze(no2_dc(:,:,41,:));
% NO2_dc_3=squeeze(no2_dc(:,:,65,:));
% NO2_dc_4=squeeze(no2_dc(:,:,t_index,:));
% 
% % Calculate the mean of the state
% 
% NO2_dc_0_mean=mean(squeeze(no2_dc(:,:,1,:)),3);
% NO2_dc_1_mean=mean(squeeze(no2_dc(:,:,17,:)),3);
% NO2_dc_2_mean=mean(squeeze(no2_dc(:,:,41,:)),3);
% NO2_dc_3_mean=mean(squeeze(no2_dc(:,:,65,:)),3);
% NO2_dc_4_mean=mean(squeeze(no2_dc(:,:,t_index,:)),3);
% 
% % Calculate the mean of the dc
% 
% NNO2_dc_0_mean=reshape(NO2_dc_0_mean,58*63,1);
% NNO2_dc_1_mean=reshape(NO2_dc_1_mean,58*63,1);
% NNO2_dc_2_mean=reshape(NO2_dc_2_mean,58*63,1);
% NNO2_dc_3_mean=reshape(NO2_dc_3_mean,58*63,1);
% NNO2_dc_4_mean=reshape(NO2_dc_4_mean,58*63,1);
% 
% 
% NNO2_dc_0=reshape(NO2_dc_0,58*63,40);
% NNO2_dc_1=reshape(NO2_dc_1,58*63,40);
% NNO2_dc_2=reshape(NO2_dc_2,58*63,40);
% NNO2_dc_3=reshape(NO2_dc_3,58*63,40);
% NNO2_dc_4=reshape(NO2_dc_4,58*63,40);
% 
% 
% NO2_state_0=squeeze(NO2_state(:,:,1,:));
% NO2_state_1=squeeze(NO2_state(:,:,17,:));
% NO2_state_2=squeeze(NO2_state(:,:,41,:));
% NO2_state_3=squeeze(NO2_state(:,:,65,:));
% NO2_state_4=squeeze(NO2_state(:,:,t_index,:));
% 
% 
% NO2_state_0_mean=mean(squeeze(NO2_state(:,:,1,:)),3);
% NO2_state_1_mean=mean(squeeze(NO2_state(:,:,17,:)),3);
% NO2_state_2_mean=mean(squeeze(NO2_state(:,:,41,:)),3);
% NO2_state_3_mean=mean(squeeze(NO2_state(:,:,65,:)),3);
% NO2_state_4_mean=mean(squeeze(NO2_state(:,:,t_index,:)),3);
%  
% 
% NNO2_state_0_mean=reshape(NO2_state_0_mean,58*63,1);
% NNO2_state_1_mean=reshape(NO2_state_1_mean,58*63,1);
% NNO2_state_2_mean=reshape(NO2_state_2_mean,58*63,1);
% NNO2_state_3_mean=reshape(NO2_state_3_mean,58*63,1);
% NNO2_state_4_mean=reshape(NO2_state_4_mean,58*63,1);
%  
% 
% 
% NNO2_state_0=reshape(NO2_state_0,58*63,40);
% NNO2_state_1=reshape(NO2_state_1,58*63,40);
% NNO2_state_2=reshape(NO2_state_2,58*63,40);
% NNO2_state_3=reshape(NO2_state_3,58*63,40);
% NNO2_state_4=reshape(NO2_state_4,58*63,40);
%  
% 
% % Calculate deviation matrix:
% for i=1:40
% 
% N_state_0_dev(:,i)=NNO2_state_0(:,i)-NNO2_state_0_mean;
% N_state_1_dev(:,i)=NNO2_state_1(:,i)-NNO2_state_1_mean;
% N_state_2_dev(:,i)=NNO2_state_2(:,i)-NNO2_state_2_mean;
% N_state_3_dev(:,i)=NNO2_state_3(:,i)-NNO2_state_3_mean;
% N_state_4_dev(:,i)=NNO2_state_4(:,i)-NNO2_state_4_mean;
% 
% N_dc_0_dev(:,i)=NNO2_dc_0(:,i)-NNO2_dc_0_mean;
% N_dc_1_dev(:,i)=NNO2_dc_1(:,i)-NNO2_dc_1_mean;
% N_dc_2_dev(:,i)=NNO2_dc_2(:,i)-NNO2_dc_2_mean;
% N_dc_3_dev(:,i)=NNO2_dc_3(:,i)-NNO2_dc_3_mean;
% N_dc_4_dev(:,i)=NNO2_dc_4(:,i)-NNO2_dc_4_mean;
% 
% end
% 
% 
% cov_0_state=(1/39)*N_state_0_dev*N_state_0_dev';
% cov_1_state=(1/39)*N_state_1_dev*N_state_1_dev';
% cov_2_state=(1/39)*N_state_2_dev*N_state_2_dev';
% cov_3_state=(1/39)*N_state_3_dev*N_state_3_dev';
% cov_4_state=(1/39)*N_state_4_dev*N_state_4_dev';
% 
% cov_0_dc=(1/39)*N_dc_0_dev*N_dc_0_dev';
% cov_1_dc=(1/39)*N_dc_1_dev*N_dc_1_dev';
% cov_2_dc=(1/39)*N_dc_2_dev*N_dc_2_dev';
% cov_3_dc=(1/39)*N_dc_3_dev*N_dc_3_dev';
% cov_4_dc=(1/39)*N_dc_4_dev*N_dc_4_dev';
% 
% 
% figure(20)
% subplot(1,3,1)
% imagesc(cov_4_state);colorbar;caxis([0 1e-19]);colormap(jet);title('L 4 Covariance EnKF standard');xlabel('states');ylabel('states')


%% Plot the state selected on the ensemble covariance matrix EnKF 
% 
% number_state_=[200*1 200*2 200*3 200*4 200*5 200*6 200*7 200*8 240*8];
% 
% figure
% subplot(2,3,1);
% 
% imagesc(cov_0_state);colorbar;caxis([0 1e-20]);colormap(jet)
% subplot(2,3,2);
% 
% imagesc(cov_1_state);colorbar;caxis([0 1e-20]);colormap(jet)
% subplot(2,3,3);
% 
% imagesc(cov_2_state);colorbar;caxis([0 1e-20]);colormap(jet)
% subplot(2,3,4);
% 
% imagesc(cov_3_state);colorbar;caxis([0 1e-20]);colormap(jet)
% subplot(2,3,5);
% 
% imagesc(cov_4_state);colorbar;caxis([0 1e-20]);colormap(jet)
% 
%% Plot the state selected on the ensemble covariance matrix Cholesky form
%  figure
%
%  imagesc(long,lati,reshape(cov_4_state(500,:),[59,63])');colorbar;caxis([0 1e-19]);colormap(cool)
% 
% 
% set(gca,'YDir','normal'); hold on
% 
% S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
% hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
% S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
% mapshow(S1,'facealpha',0)
%   
% % S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% % mapshow(S2,'facealpha',0)
% 
% xlim([-76 -71.9]);ylim([9 13.5]);h=colorbar;ylabel(h,'1e15 mlc/cm2');title('Covariance EnKF state 500');xlabel('longitude °');ylabel('latitude °')
% % saveas(fig,'Tropomi Retrieval','jpg')
% 
% %
% colorbar;
% 
% 
% figure
% for i=1:9
% number_state=number_state_(i);lati_state=number_state/59;
% subplot(3,3,i)
% imagesc(long,lati,reshape(cov_4_state(number_state,:),[59,63])');colorbar;caxis([0 1e-19]);colormap(cool)
% hold on
% s=scatter(long(number_state-59*floor(lati_state)),lati(floor(lati_state)),100,'s','MarkerEdgeColor',[1 1 1],...
%               'MarkerFaceColor',[ 0 0 0],...
%               'LineWidth',1.5,'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8); 
% set(gca,'YDir','normal'); hold on
% S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
% hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
% S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
% mapshow(S1,'facealpha',0)
% % S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% % mapshow(S2,'facealpha',0)
% xlim([-76 -71.9]);ylim([9 13.5]);h=colorbar;ylabel(h,'1e15 mlc/cm2');title(sprintf('Covariance EnKF state= %i',number_state_(i)));xlabel('longitude °');ylabel('latitude °')
% % saveas(fig,'Tropomi Retrieval','jpg')
% colorbar;
% end
% % 
% %% Concatenate states + dc
% 
% figure
% state_dc_lag0=[NNO2_state_0;NNO2_dc_0(1,:)];
% state_dc_lag1=[NNO2_state_1;NNO2_dc_1(1,:)];
% state_dc_lag2=[NNO2_state_2;NNO2_dc_2(1,:)];
% state_dc_lag3=[NNO2_state_3;NNO2_dc_3(1,:)];
% state_dc_lag4=[NNO2_state_4;NNO2_dc_4(1,:)];
% 
% state_lag4_state=[NNO2_state_4];
% 
% % B_state_dc_lag0=Calculo_B_Cholesky(state_dc_lag0,1);
% % B_state_dc_lag1=Calculo_B_Cholesky(state_dc_lag1,1);
% % B_state_dc_lag2=Calculo_B_Cholesky(state_dc_lag2,1);
% % B_state_dc_lag3=Calculo_B_Cholesky(state_dc_lag3,1);
%  B_state_dc_lag4=Calculo_B_Cholesky(state_dc_lag4,1);
%  B_state_state_lag4=Calculo_B_Cholesky(state_lag4_state,1);
%  
% 
% % subplot(2,3,1);imagesc(B_state_dc_lag0);colorbar;caxis([0 1e-15]);colormap(jet);title('L 0');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% % subplot(2,3,2);imagesc(B_state_dc_lag1);colorbar;caxis([0 1e-22]);colormap(jet);title('L 1');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% % subplot(2,3,3);imagesc(B_state_dc_lag2);colorbar;caxis([0 1e-23]);colormap(jet);title('L 2');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% % subplot(2,3,4);imagesc(B_state_dc_lag3);colorbar;caxis([0 1e-24]);colormap(jet);title('L 3');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% % subplot(2,3,5);imagesc(B_state_dc_lag4);colorbar;caxis([0 1e-25]);colormap(jet);title('L 4');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% % 
% 
%  
% figure
% subplot(1,2,1)
% imagesc(normalize_beta(B_state_dc_lag4));colorbar;%caxis([0 1e-16]);colormap(jet);title('L 4 Covariance Modified Cholesky ');xlabel('states');ylabel('states')
% subplot(1,2,2)
% imagesc(normalize_beta(B_state_state_lag4));colorbar;%caxis([0 1e-16]);colormap(jet);title('L 4 Covariance Modified Cholesky ');xlabel('states');ylabel('states')
% 
%  
%% 

% number_state_=[200 201 202 203 204 205 206 207 300];
% 
% 
% figure
% for i=1:9
% number_state=number_state_(i);lati_state=floor(number_state/58);
% subplot(3,3,i)
% imagesc(long,lati,reshape(B_state_dc_lag4(number_state,:),[58,63])');colorbar;caxis([0 1e-17]);colormap(cool)
% hold on
% s=scatter(long(number_state-58*floor(lati_state)),lati(floor(lati_state)),100,'s','MarkerEdgeColor',[1 1 1],...
%               'MarkerFaceColor',[ 0 0 0],...
%               'LineWidth',1.5,'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8); 
% set(gca,'YDir','normal'); hold on
% S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
% hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
% S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
% mapshow(S1,'facealpha',0)
% % S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% % mapshow(S2,'facealpha',0)
% xlim([-76 -71.9]);ylim([9 13.5]);h=colorbar;ylabel(h,'1e15 mlc/cm2');title(sprintf('Covariance MC state= %i',number_state_(i)));xlabel('longitude °');ylabel('latitude °')
% % saveas(fig,'Tropomi Retrieval','jpg')
% colorbar;
% end

%%  Calculate Analysys step

% X_a=X_b+sum(LD^(1/2)Y^T[YY^T]d_i


% X_a=X_b+sum((B^(-1)+H^(T)RH)^1)HRd_i)


%%  Part to calculate Y=HX~ [Y_1-Ym,Y_2-Ym,...,Y_N-Ym]   siendo N el número de ensambles
figure

% move to the folder where the outputs are

mydir='/run/media/dirac/Datos/scratch/projects/Prueba_numero_4_EnKS_MC/Prueba_numero_4_EnKS_MC/output';
cd(mydir)

% name of the outputs of interest

files_y_ensembles={'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190201_1800_xi','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190202_1800_xi',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1800_xi','LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190203_1900_xi',...
    'LE_Prueba_numero_4_EnKS_MC_tropomi-no2_20190204_1900_xi'};


% Import values in a matrix remebering to substract the validation value

ens=40;
tic
for j=1:5   % each day of the window collecting the ensemble observation output
subplot(1,5,j)
ilat=ncread(files{j},'ilat');ilon=ncread(files{j},'ilon');
  V=cell2mat(randam_obs(j)); 
  
for i=1:ens
  
    if(i<10)
        y_ens=ncread(strcat(files_y_ensembles{j},sprintf('0%ia.nc',i)),'y');
        Validate_ENS(i,:)=y_ens(V);
        y_ens(V)=[];
        Y_ENS_(i,:)=y_ens;   
    end
    
     if(i>=10)
         y_ens=ncread(strcat(files_y_ensembles{j},sprintf('%ia.nc',i)),'y');
         Validate_ENS(i,:)=y_ens(V);
         y_ens(V)=[];
         Y_ENS_(i,:)=y_ens;
     end

    
     plot(y_ens);hold on

end
Almacenando_asim{j,:,:}=Y_ENS_;
Almacenando_Validate{j,:,:}=Validate_ENS;
clear V Y_ENS_ y_ens Validate_ENS
end
% calculate the mean
time=toc
% calculate


save('/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/Y_ENS.mat','Almacenando_asim','Almacenando_Validate')







 


