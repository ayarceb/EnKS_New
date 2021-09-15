% Code to implement the EnKS_MC technique with the LOTOS-EUROS for the
% north region Colombia domain
clc;close all;clear all
addpath('/home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC_new')  
% Variable characteristics    Size no2:       58x63x1x97
%                             Dimensions: longitude,latitude,level,time
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


%no2=zeros(45,52,97,40); no2_dc=zeros(45,52,97,40);
no2=zeros(58,63,97,40); no2_dc=zeros(58,63,97,40);

t1 = datetime(2019,2,1,0,0,0);t2 = datetime(2019,2,5,0,0,0);t = t1:hours(1):t2;


ens=40; %NO2_sfc Ensemble number concatenate
for j=1:1
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

%% take state snapshots from the forward model

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




%% Different Lags and having the radious constant

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
[B_chol,B_square]=Calculo_B_Cholesky(STATES_lag4,1);
%%
% figure

% subplot(2,3,1);imagesc(B_lag0);colorbar;colormap(jet);caxis([0 1e-17]);title('L 0');xlabel('states');ylabel('states')
% subplot(2,3,2);imagesc(B_lag1);colorbar;colormap(jet);caxis([0 1e-17]);title('L 1');xlabel('states');ylabel('states')
% subplot(2,3,3);imagesc(B_lag2);colorbar;colormap(jet);caxis([0 1e-17]);title('L 2');xlabel('states');ylabel('states')
% subplot(2,3,4);imagesc(B_lag3);colorbar;colormap(jet);caxis([0 1e-17]);title('L 3');xlabel('states');ylabel('states')

% subplot(2,3,5);
% imagesc(B_lag4);colorbar;colormap(jet);caxis([0 1e-17]);title('L 4');xlabel('states');ylabel('states')
%%  Plot Cholesky and Square Cholesky
figure
subplot(1,2,1) 
imagesc(B_square)
colorbar;colormap(jet)
caxis([0 1e-9]);ylabel('States');xlabel('States')
title('Square of the Cholesky matrix')
subplot(1,2,2)
imagesc(B_chol)
colorbar;colormap(jet)
caxis([0 1e-16]);ylabel('States');xlabel('States')
title('Complete Cholesky matrix')


%% PLOT DC

figure

for j=1:1
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

%% take snapshots of the DC factor from the forward model

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
% parameters
% t_index=97;
% 
% NO2_dc_0=squeeze(no2_dc(:,:,1,:));
% NO2_dc_1=squeeze(no2_dc(:,:,17,:));
% NO2_dc_2=squeeze(no2_dc(:,:,41,:));
% NO2_dc_3=squeeze(no2_dc(:,:,65,:));
% NO2_dc_4=squeeze(no2_dc(:,:,t_index,:));
% 
% NO2_dc_0_mean=mean(squeeze(no2_dc(:,:,1,:)),3);
% NO2_dc_1_mean=mean(squeeze(no2_dc(:,:,17,:)),3);
% NO2_dc_2_mean=mean(squeeze(no2_dc(:,:,41,:)),3);
% NO2_dc_3_mean=mean(squeeze(no2_dc(:,:,65,:)),3);
% NO2_dc_4_mean=mean(squeeze(no2_dc(:,:,t_index,:)),3);
% 
% 
% NNO2_dc_0_mean=reshape(NO2_dc_0_mean,58*63,1);
% NNO2_dc_1_mean=reshape(NO2_dc_1_mean,58*63,1);
% NNO2_dc_2_mean=reshape(NO2_dc_2_mean,58*63,1);
% NNO2_dc_3_mean=reshape(NO2_dc_3_mean,58*63,1);
% NNO2_dc_4_mean=reshape(NO2_dc_4_mean,58*63,1);
% 
% 
% 
% NNO2_dc_0=reshape(NO2_dc_0,58*63,40);
% NNO2_dc_1=reshape(NO2_dc_1,58*63,40);
% NNO2_dc_2=reshape(NO2_dc_2,58*63,40);
% NNO2_dc_3=reshape(NO2_dc_3,58*63,40);
% NNO2_dc_4=reshape(NO2_dc_4,58*63,40);
% 
% %states
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
% 
% NNO2_state_0_mean=reshape(NO2_state_0_mean,58*63,1);
% NNO2_state_1_mean=reshape(NO2_state_1_mean,58*63,1);
% NNO2_state_2_mean=reshape(NO2_state_2_mean,58*63,1);
% NNO2_state_3_mean=reshape(NO2_state_3_mean,58*63,1);
% NNO2_state_4_mean=reshape(NO2_state_4_mean,58*63,1);
% 
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
% 
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
% %% Plot the state selected on the ensemble covariance matrix EnKF 
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
% 
% figure
% 
% imagesc(cov_4_state);colorbar;caxis([0 1e-19]);colormap(jet);title('L 4 Covariance EnKF standard');xlabel('states');ylabel('states')
% 
% 
% figure
% %     i
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


%%

% figure
% state_dc_lag0=[NNO2_state_0;NNO2_dc_0(1,:)];
% state_dc_lag1=[NNO2_state_0;NNO2_dc_0(1,:);NNO2_state_1;NNO2_dc_1(1,:)];
% state_dc_lag2=[NNO2_state_0;NNO2_dc_0(1,:);NNO2_state_1;NNO2_dc_1(1,:);NNO2_state_2;NNO2_dc_2(1,:)];
% state_dc_lag3=[NNO2_state_0;NNO2_dc_0(1,:);NNO2_state_1;NNO2_dc_1(1,:);NNO2_state_2;NNO2_dc_2(1,:);NNO2_state_3;NNO2_dc_3(1,:)];
% state_dc_lag4=[NNO2_state_0;NNO2_dc_0(1,:);NNO2_state_1;NNO2_dc_1(1,:);NNO2_state_2;NNO2_dc_2(1,:);NNO2_state_3;NNO2_dc_3(1,:);NNO2_state_4;NNO2_dc_4(1,:)];
% 
% 
%  B_state_dc_lag0=Calculo_B_Cholesky(state_dc_lag0,1);
%  B_state_dc_lag1=Calculo_B_Cholesky(state_dc_lag1,1);
%  B_state_dc_lag2=Calculo_B_Cholesky(state_dc_lag2,1);
%  B_state_dc_lag3=Calculo_B_Cholesky(state_dc_lag3,1);
%  B_state_dc_lag4=Calculo_B_Cholesky(state_dc_lag4,1);

% subplot(2,3,1);imagesc(B_state_dc_lag0);colorbar;caxis([0 1e-15]);colormap(jet);title('L 0');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% subplot(2,3,2);imagesc(B_state_dc_lag1);colorbar;caxis([0 1e-22]);colormap(jet);title('L 1');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% subplot(2,3,3);imagesc(B_state_dc_lag2);colorbar;caxis([0 1e-23]);colormap(jet);title('L 2');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% subplot(2,3,4);imagesc(B_state_dc_lag3);colorbar;caxis([0 1e-24]);colormap(jet);title('L 3');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')
% subplot(2,3,5);imagesc(B_state_dc_lag4);colorbar;caxis([0 1e-25]);colormap(jet);title('L 4');xlabel('states + DC Factor parameter');ylabel('states + DC Factor parameter')

% states_lag0=[NNO2_state_0];
% states_lag1=[NNO2_state_0;NNO2_state_1];
% states_lag2=[NNO2_state_0;NNO2_state_1;NNO2_state_2];
% states_lag3=[NNO2_state_0;NNO2_state_1;NNO2_state_2;NNO2_state_3];
% states_lag4=[NNO2_state_0;NNO2_state_1;NNO2_state_2;NNO2_state_3;NNO2_state_4];

states_lag0=[NNO2_state_0];
states_lag1=[NNO2_state_1];
states_lag2=[NNO2_state_2];
states_lag3=[NNO2_state_3];
states_lag4=[NNO2_state_4];



%  B_state_lag0=Calculo_B_Cholesky(states_lag0,1);
%  B_state_lag1=Calculo_B_Cholesky(states_lag1,1);
%  B_state_lag2=Calculo_B_Cholesky(states_lag2,1);
%  B_state_lag3=Calculo_B_Cholesky(states_lag3,1);
 B_state_lag4=Calculo_B_Cholesky(NNO2_state_4,1);
 
 
figure

imagesc(B_state_lag4);colorbar;caxis([0 1e-16]);colormap(jet);title('L 4 Covariance Modified Cholesky ');xlabel('states');ylabel('states')

 
%% 

number_state_=[200 201 202 203 204 205 206 207 300];


figure
for i=1:9
number_state=number_state_(i);lati_state=floor(number_state/58);
subplot(3,3,i)
imagesc(long,lati,reshape(B_state_lag4(number_state,:),[58,63])');colorbar;caxis([0 1e-17]);colormap(cool)
hold on
s=scatter(long(number_state-58*floor(lati_state)),lati(floor(lati_state)),100,'s','MarkerEdgeColor',[1 1 1],...
              'MarkerFaceColor',[ 0 0 0],...
              'LineWidth',1.5,'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8); 
set(gca,'YDir','normal'); hold on
S=shaperead('/run/media/dirac/Datos/Real_DROPBOX/Dropbox/2017/Doctorado/SIG/05_ANTIOQUIA_/ADMINISTRATIVO/MGN_ADM_DPTO_POLITICO.shp');
hold on; mapshow(S,'FaceAlpha',0, 'LineWidth',1)
S1=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/MAUI.LatinAmerica.EPSG4326/MAUI.LatinAmerica.EPSG4326.shp')
mapshow(S1,'facealpha',0)
% S2=shaperead('/run/media/dirac/Datos/Reciente_Dropbox/2018/SIG/gadm36_PAN_shp/gadm36_PAN_0.shp')
% mapshow(S2,'facealpha',0)
xlim([-76 -71.9]);ylim([9 13.5]);h=colorbar;ylabel(h,'1e15 mlc/cm2');title(sprintf('Covariance MC state= %i',number_state_(i)));xlabel('longitude °');ylabel('latitude °')
% saveas(fig,'Tropomi Retrieval','jpg')
colorbar;
end



