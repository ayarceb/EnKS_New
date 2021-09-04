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

% Coordinates for the 59x63 array.
% Coordenadas        (latitude,longitude)
% Barranquilla       Latitud: 10.9878,    Longitud: -74.7889      (33,20)
% Santa Marta        Latitud: 11.24079,   Longitud: -74.19904      (37,26)
% Cartagena          Latitud: 10.3884731, Longitud: -75.488593     (27,12)
% Mina               Latitud: 9.5829857,  Longitud: -73.479376      (18,35)
% Valledupar         Latitud: 10.4646054, Longitud: -73.2599952     ( 27,36)

ciudades={'Barranquilla ','Santa Marta ','Cartagena ','Mina Drummond ','Valledupar '};

lon=[20,26,12,35,36];lat=[33,37,27,18,27];

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

% take state snapshots from the forward model

% NO2_state_0=squeeze(NO2_state(:,:,1,:));
% NO2_state_1=squeeze(NO2_state(:,:,17,:));
% NO2_state_2=squeeze(NO2_state(:,:,41,:));
% NO2_state_3=squeeze(NO2_state(:,:,65,:));
% NO2_state_4=squeeze(NO2_state(:,:,89,:));
% 
% 
% NNO2_state_0=reshape(NO2_state_0,45*52,40);
% NNO2_state_1=reshape(NO2_state_1,45*52,40);
% NNO2_state_2=reshape(NO2_state_2,45*52,40);
% NNO2_state_3=reshape(NO2_state_3,45*52,40);
% NNO2_state_4=reshape(NO2_state_4,45*52,40);


%% Different radious with the Modified Cholesky plots
% STATES=[NNO2_state_0;NNO2_state_1;NNO2_state_2;NNO2_state_3;NNO2_state_4];
% 
% 
% 
% Cov1=cov(STATES');
 !% Cov2=cov(STATES*STATES');
% 
% B1=Calculo_B_Cholesky(STATES,1);
% B2=Calculo_B_Cholesky(STATES,2);
% B3=Calculo_B_Cholesky(STATES,3);
% B4=Calculo_B_Cholesky(STATES,4);

%% Different Lags and having the radious constant

% figure
% STATES_lag0=[NNO2_state_0];
% STATES_lag1=[NNO2_state_0;NNO2_state_1];
% STATES_lag2=[NNO2_state_0;NNO2_state_1;NNO2_state_2];
% STATES_lag3=[NNO2_state_0;NNO2_state_1;NNO2_state_2;NNO2_state_3];
% STATES_lag4=[NNO2_state_0;NNO2_state_1;NNO2_state_2;NNO2_state_3;NNO2_state_4];
% 
% 
% B_lag0=Calculo_B_Cholesky(STATES_lag0,1);
% B_lag1=Calculo_B_Cholesky(STATES_lag1,1);
% B_lag2=Calculo_B_Cholesky(STATES_lag2,1);
% B_lag3=Calculo_B_Cholesky(STATES_lag3,1);
% B_lag4=Calculo_B_Cholesky(STATES_lag4,1);

%%
% figure
% 
% subplot(2,3,1);imagesc(B_lag0);colorbar;colormap(jet);caxis([0 1e-17]);title('L 0');xlabel('states');ylabel('states')
% subplot(2,3,2);imagesc(B_lag1);colorbar;colormap(jet);caxis([0 1e-17]);title('L 1');xlabel('states');ylabel('states')
% subplot(2,3,3);imagesc(B_lag2);colorbar;colormap(jet);caxis([0 1e-17]);title('L 2');xlabel('states');ylabel('states')
% subplot(2,3,4);imagesc(B_lag3);colorbar;colormap(jet);caxis([0 1e-17]);title('L 3');xlabel('states');ylabel('states')
% subplot(2,3,5);imagesc(B_lag4);colorbar;colormap(jet);caxis([0 1e-17]);title('L 4');xlabel('states');ylabel('states')





%% Plots Covariance vs Cholesky


% figure
% 
% subplot(2,3,1);imagesc(Cov1);colorbar;colormap(jet);caxis([0 1e-19]);title('cov(STATES)');xlabel('states');ylabel('states')
% subplot(2,3,2);imagesc(Cov2);colorbar;colormap(jet);caxis([0 1e-33]);title('cov(STATES)');xlabel('states');ylabel('states')
% subplot(2,3,3);imagesc(B1);colorbar;colormap(jet);caxis([0 1e-17]);title('Cholesky r=1');xlabel('states');ylabel('states')
% subplot(2,3,4);imagesc(B2);colorbar;colormap(jet);caxis([0 1e-16]);title('Cholesky r=2');xlabel('states');ylabel('states')
% subplot(2,3,5);imagesc(B3);colorbar;colormap(jet);caxis([0 1e-15]);title('Cholesky r=3');xlabel('states');ylabel('states')
% subplot(2,3,6);imagesc(B4);colorbar;colormap(jet);caxis([0 1e-14]);title('Cholesky r=4');xlabel('states');ylabel('states')


%% With the parameters




%% PLOT DC

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
mydir=append('/run/media/dirac/Datos/scratch/projects/',name_run,'/',name_run,'/output');cd(mydir)

% %% take snapshots of the DC factor from the forward model
% figure
% NO2_dc_0=squeeze(no2_dc(:,:,1,:));
% NO2_dc_1=squeeze(no2_dc(:,:,17,:));
% NO2_dc_2=squeeze(no2_dc(:,:,41,:));
% NO2_dc_3=squeeze(no2_dc(:,:,65,:));
% NO2_dc_4=squeeze(no2_dc(:,:,89,:));
% 
% 
% NNO2_dc_0=reshape(NO2_dc_0,45*52,40);
% NNO2_dc_1=reshape(NO2_dc_1,45*52,40);
% NNO2_dc_2=reshape(NO2_dc_2,45*52,40);
% NNO2_dc_3=reshape(NO2_dc_3,45*52,40);
% NNO2_dc_4=reshape(NO2_dc_4,45*52,40);
% 
% figure
% dc_lag0=[NNO2_dc_0(1,:)];
% dc_lag1=[NNO2_dc_0(1,:);NNO2_dc_1(1,:)];
% dc_lag2=[NNO2_dc_0(1,:);NNO2_dc_1(1,:);NNO2_dc_2(1,:)];
% dc_lag3=[NNO2_dc_0(1,:);NNO2_dc_1(1,:);NNO2_dc_2(1,:);NNO2_dc_3(1,:)];
% dc_lag4=[NNO2_dc_0(1,:);NNO2_dc_1(1,:);NNO2_dc_2(1,:);NNO2_dc_3(1,:);NNO2_dc_4(1,:)];
% 
% B_dc_lag0=Calculo_B_Cholesky(dc_lag0,1);
% B_dc_lag1=Calculo_B_Cholesky(dc_lag1,1);
% B_dc_lag2=Calculo_B_Cholesky(dc_lag2,1);
% B_dc_lag3=Calculo_B_Cholesky(dc_lag3,1);
% B_dc_lag4=Calculo_B_Cholesky(dc_lag4,1);
% 
% subplot(2,3,1);imagesc(B_dc_lag0);colorbar;colormap(jet);title('L 0');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% subplot(2,3,2);imagesc(B_dc_lag1);colorbar;colormap(jet);title('L 1');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% subplot(2,3,3);imagesc(B_dc_lag2);colorbar;colormap(jet);title('L 2');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% subplot(2,3,4);imagesc(B_dc_lag3);colorbar;colormap(jet);title('L 3');xlabel('DC Factor parameter');ylabel('DC Factor parameter')
% subplot(2,3,5);imagesc(B_dc_lag4);colorbar;colormap(jet);title('L 4');xlabel('DC Factor parameter');ylabel('DC Factor parameter')

%% States + parameteres

NO2_dc_0=squeeze(no2_dc(:,:,1,:));
NO2_dc_1=squeeze(no2_dc(:,:,17,:));
NO2_dc_2=squeeze(no2_dc(:,:,41,:));
NO2_dc_3=squeeze(no2_dc(:,:,65,:));
NO2_dc_4=squeeze(no2_dc(:,:,89,:));

NNO2_dc_0=reshape(NO2_dc_0,45*52,40);
NNO2_dc_1=reshape(NO2_dc_1,45*52,40);
NNO2_dc_2=reshape(NO2_dc_2,45*52,40);
NNO2_dc_3=reshape(NO2_dc_3,45*52,40);
NNO2_dc_4=reshape(NO2_dc_4,45*52,40);


NO2_state_0=squeeze(NO2_state(:,:,1,:));
NO2_state_1=squeeze(NO2_state(:,:,17,:));
NO2_state_2=squeeze(NO2_state(:,:,41,:));
NO2_state_3=squeeze(NO2_state(:,:,65,:));
NO2_state_4=squeeze(NO2_state(:,:,89,:));

NNO2_state_0=reshape(NO2_state_0,45*52,40);
NNO2_state_1=reshape(NO2_state_1,45*52,40);
NNO2_state_2=reshape(NO2_state_2,45*52,40);
NNO2_state_3=reshape(NO2_state_3,45*52,40);
NNO2_state_4=reshape(NO2_state_4,45*52,40);



figure
state_dc_lag0=[NNO2_state_0;NNO2_dc_0(1,:)];
state_dc_lag1=[NNO2_state_0;NNO2_dc_0(1,:);NNO2_state_1;NNO2_dc_1(1,:)];
state_dc_lag2=[NNO2_state_0;NNO2_dc_0(1,:);NNO2_state_1;NNO2_dc_1(1,:);NNO2_state_2;NNO2_dc_2(1,:)];
state_dc_lag3=[NNO2_state_0;NNO2_dc_0(1,:);NNO2_state_1;NNO2_dc_1(1,:);NNO2_state_2;NNO2_dc_2(1,:);NNO2_state_3;NNO2_dc_3(1,:)];
state_dc_lag4=[NNO2_state_0;NNO2_dc_0(1,:);NNO2_state_1;NNO2_dc_1(1,:);NNO2_state_2;NNO2_dc_2(1,:);NNO2_state_3;NNO2_dc_3(1,:);NNO2_state_4;NNO2_dc_4(1,:)];


% B_state_dc_lag0=Calculo_B_Cholesky(state_dc_lag0,1);
% B_state_dc_lag1=Calculo_B_Cholesky(state_dc_lag1,1);
% B_state_dc_lag2=Calculo_B_Cholesky(state_dc_lag2,1);
%B_state_dc_lag3=Calculo_B_Cholesky(state_dc_lag3,1);
B_state_dc_lag4=Calculo_B_Cholesky(state_dc_lag4,1);


subplot(2,3,1);imagesc(B_state_dc_lag0);colorbar;caxis([0 1e-15]);colormap(jet);title('L 0');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
subplot(2,3,2);imagesc(B_state_dc_lag1);colorbar;caxis([0 1e-22]);colormap(jet);title('L 1');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
subplot(2,3,3);imagesc(B_state_dc_lag2);colorbar;caxis([0 1e-23]);colormap(jet);title('L 2');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
subplot(2,3,4);imagesc(B_state_dc_lag3);colorbar;caxis([0 1e-24]);colormap(jet);title('L 3');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
subplot(2,3,5);imagesc(B_state_dc_lag4);colorbar;caxis([0 1e-25]);colormap(jet);title('L 4');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')


%%

figure(1)
Xb=nan(45,52,97,40);
Fe=nan(45,52,97,40);

Xb_mean=nan(5,97);
Fe_mean=nan(5,97);

%% This segment put the points for the time snapshot selected which correspond to the TROPOMI time

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
% figure(2)
% for j=1:5
% subplot(2,3,j)
%    hhh= plot(t(1:end-4),squeeze(Fe(lon(j),lat(j),5:end,1)),'Marker', 'o', 'MarkerSize', 4, ...
%       'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%   ppp= plot(t(1:end-4),squeeze(Fe_mean(j,5:end)),'Marker', 'o', 'MarkerSize', 4, ...
%      'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%  legend([h,p,hhh,ppp],'Ensemble members','Mean','Background ensemble member state', 'Background mean state')
% end