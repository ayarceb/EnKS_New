% Code to run the experiment for the LOTOS-EUROS EnKS-MC
% This code run from matlab the LOTOS-EUROS performing a EnKS-MC Data
% assimilation scheme
% first task: Assimilate with the synthetic NO2 columns


% the bash file RunLEKF.sh configurate the run and launch the lekf with the emission factor perturbation 


clc;close all;clear all


%%===Path where the sh file to run the ensembles is located===
mydir='/home/dirac/Dropbox/2020/EnKS-MC/EnKS-MC/ENKS_MC_LOTOS'
cd(mydir)
!./RunLEKF.sh



%% Merge outputs
cd()



%% 


     meanxb(:,i)=mean(X_S(:,:,i),2);    
    
  for Nen=1:N
     XXb(:,Nen,i)= X_S(:,Nen,i)-meanxb(:,i);  
  end
  
%  Standard EnKF Background matrix 

   B(:,:,i)=((1/N-1)*(XXb(:,:,i)*XXb(:,:,i)'));








%% %==========Scenario EnKS ==================
% if select==1
%   
for i=2:Tsim-M*frequency
%======== Forecast Step==========

    for en=1:N
         [X_S(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(Xa(:,en,i-1)),F);
    end
  
     meanxb(:,i)=mean(X_S(:,:,i),2);    
    
  for Nen=1:N
     XXb(:,Nen,i)= X_S(:,Nen,i)-meanxb(:,i);  
  end
  
%  Standard EnKF Background matrix 

   B(:,:,i)=((1/N-1)*(XXb(:,:,i)*XXb(:,:,i)'));

   if sum(muestreo==i)   %Validation if there are a observation available
%======================== Analysis Step ENKS=========================================

    X_aw(:,:,1)=X_S(:,:,i);
    
    if M==1
    B_aw(:,:,1)=B(:,:,i);    
    end
        for kk=2:(M*frequency)
            for en=1:N    
                [X_aw(:,en,kk)]=Lorenz_96_one_step(1,dt,squeeze(X_aw(:,en,kk-1)),F);
            end
            meanxb_aw(:,kk)=mean(X_aw(:,:,kk),2);    
    
             for Nen=1:N
                XXb_aw(:,Nen,kk)= X_aw(:,Nen,kk)-meanxb_aw(:,kk)+1*randn(n,1);  
             end
  
             B_aw(:,:,kk)=((1/N-1)*(XXb_aw(:,:,kk)*XXb_aw(:,:,kk)'));   
        end
    aux=0;
    for kk2=1:frequency:M*frequency
        K(:,:,kk2)=B_aw(:,:,kk2)*H'*pinv(H*B_aw(:,:,kk2)*H'+R);     % Matriz de Ganancia de Kalman version 
        C(:,kk2)= Y(:,i+kk2-1)+sigma*randn(m,1);
        aux=aux+(K(:,:,kk2)*(C(:,kk2)-H*X_aw(:,:,kk2)))/M;
    end
    
    for en=1:N
      Xa(:,en,i)=X_S(:,en,i)+aux(:,en);
    end
    

else 
  Xa(:,:,i)=X_S(:,:,i);  
end
  meanxa_EnKS(:,i)=mean(Xa(:,:,i),2);                                                                                                                                                                                                                                                                                                                                              
end

%% ================================Scenario EnKS-MC ============================================================================================================================
% %% 
% % if select==1
% %   
 for i=2:Tsim-M*frequency
% %======== Forecast Step==========
 i
     for en=1:N
          [X_S_chol(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(X_S_chol(:,en,i-1)),F);
     end
   
      meanxb(:,i)=mean(X_S_chol(:,:,i),2);    
     
   for Nen=1:N
      XXb(:,Nen,i)= X_S_chol(:,Nen,i)-meanxb(:,i);  
   end
   
 %  Standard Modified cholesky Background matrix 
% 
     % BC(:,:,i)=Calculo_B_Cholesky(XXb,r); %Estimation of Covariance by Modified Cholesky
      BC(:,:,i)=pinv(Calculo_B_Cholesky(XXb,r)); %Estimation of Covariance by Modified Cholesky
      

 if sum(muestreo==i)   %Validation if there are a observation available
% %======================== Analysis Step ENKS=========================================
% 
     X_aw(:,:,1)=X_S_chol(:,:,i);
     
     if M==1
     B_aw(:,:,1)=BC(:,:,i);    
     end
         for kk=2:(M*frequency)
             for en=1:N    
                 [X_aw(:,en,kk)]=Lorenz_96_one_step(1,dt,squeeze(X_aw(:,en,kk-1)),F);
             end
             meanxb_aw(:,kk)=mean(X_aw(:,:,kk),2);    
     
              for Nen=1:N
                 XXb_aw(:,Nen,kk)= X_aw(:,Nen,kk)-meanxb_aw(:,kk)+0.2*randn(n,1);  
              end
   
               B_aw(:,:,kk)=pinv(Calculo_B_Cholesky(XXb_aw,r));   
         end
     aux=0;
     for kk2=1:frequency:M*frequency
         K(:,:,kk2)=B_aw(:,:,kk2)*H'*pinv(H*B_aw(:,:,kk2)*H'+R);     % Matriz de Ganancia de Kalman
         C(:,kk2)= Y(:,i+kk2-1)+sigma*randn(m,1);
         aux=aux+(K(:,:,kk2)*(C(:,kk2)-H*X_aw(:,:,kk2)))/M;
     end
%     
     for en=1:N
       X_S_chol(:,en,i)=X_S_chol(:,en,i)+aux(:,en);
     end
     
% 
 else 
   X_S_chol(:,:,i)=X_S_chol(:,:,i);  
 end
   meanxa_EnKS_chol(:,i)=mean(X_S_chol(:,:,i),2);                                                                                                                                                                                                                                                                                                                                              
 end
% 