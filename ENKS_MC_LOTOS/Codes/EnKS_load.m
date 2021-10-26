%  Code to load and calculate the analysis step of the EnKS_MC
%  Andrés Yarce Botero 2021-02

%  The matrices loaded to this code comes from the code:
%  /home/dirac/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC_new/ENKS_MC_LOTOS/Codes/EnKS_preparing_V2.m

%  Calculate Analysys step

% X_a=X_b+sum(LD^(1/2)(Y^T[YY^T])^(-1)d_i
% X_a=X_b+sum((B^(-1)+H^(T)RH)^1)HRd_i)


clc;close all;clear all

load /run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/Matrices.mat
load /run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/innovation.mat
load /run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/Y_ENS.mat
load /run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/Sigmas.mat
load /run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/ENKS_MC/save_var/B_square.mat


for j=1:5
mean_Y_Ens_assim{j,:}=mean(Almacenando_asim{j}); 

for i=1:40
      Aux=cell2mat(Almacenando_asim(j));
    Y(i,:)=Aux(i,:)-cell2mat(mean_Y_Ens_assim(j));
    
end

inno=cell2mat(d(j));
Y*[Y'*Y]*inno'

Y_aux{j,:,:}=Y;

 clear Aux Y
end


