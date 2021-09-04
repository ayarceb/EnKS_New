% EnKF script for the Lorenz96 using th Modified Cholesky approximation of
% the B covariance matrix

clear all
close all
clc

%%

Tsim=200;        % Simulation time
dt=0.01;         % Step length
m=40;            % Number of observations
n=40;            % State number
F=8;             % Forcing factor model
r=10;             % Predecesors radio
sigma=0.5;
R=sigma^2*eye(m);       % Observation Covariance matrix 
H = eye(n,n);           % Observation operator matrix
H = H(randperm(n,m),:);

%===Generate real state===
x0=1*rand(n,1);
[Xreal]=Lorenz_96(Tsim,dt,x0,F);            % Create vecto


%==Number of Ensembles==
N=100;                    % Number of ensembles
Xb=zeros(n,N,Tsim);       % Background ensemble initialization
meanxb=zeros(n,Tsim);       % Background ensemble initialization
Xb(:,:,1)=1*rand(n,N);    
Xa=Xb;

meanxa_EnKF=zeros(n,Tsim);
meanxa_EnKF_CHOLESKY=zeros(n,Tsim);


% %%
% %=====Scenario EnKF===
% %==Observations===
% Y=H*Xreal;     %  Sampling of the observation operator
% for i=1:Tsim-1
%   %===== Forecast Step=====
%   for en=1:N
%      [Xb(:,en,i+1)]=Lorenz_96_one_step(1,dt,squeeze(Xb(:,en,i)),F);
%   end
%      meanxb(:,i)=mean(Xb(:,:,i),2);    
%   for Nen=1:N
%      XXb(:,Nen,i)= Xb(:,Nen,i)-meanxb(:,i);  
%   end
%   i
% 
%    %XXa=XXb(:,:,i)*sqrt(1/(N-1));   
%    BC=Calculo_B_Cholesky(XXb,r); %Estimation of Covariance by Modified Cholesky
% %   inB=Calculo_inB_Cholesky(XXa,r); %Estimation of invserse Covariance by Modified Cholesky
% %   Bsquare=B^(1/2);
% 
%    
%      
% %     L=Xb(:,:,i)-meanxb;
% %     
% %%  Standard EnKF Background matrix 
% 
%    B=((1/N-1)*(XXb(:,:,i)*XXb(:,:,i)'));
% 
% 
% %============== Analysis Step ENKF=========================================
% muestra=0;     
% K=B*H'*pinv(H*B*H'+R);     % Matriz de Ganancia de Kalman
%      
%      
%      for en=1:N
%          
%          C(:,i)= Y(:,i)+sigma*randn(m,1);        %  Almacenamiento de los datos sintéticos
%          Xa(:,en,i+1)=Xb(:,en,i+1)+K*(C(:,i)-H*Xb(:,en,i+1));
% 
%      end
%      meanxa_EnKF(:,i+1)=mean(Xa(:,:,i+1),2);
%      
% %=============== Analysis Step ENKF========================================     
%      KC=BC*H'*pinv(H*BC*H'+R);     % Matriz de Ganancia de Kalman para cholesky
%      
%      
%      for en=1:N
%          
%          C(:,i)= Y(:,i)+sigma*randn(m,1);        %  Almacenamiento de los datos sintéticos
%          Xa_C(:,en,i+1)=Xb(:,en,i+1)+KC*(C(:,i)-H*Xb(:,en,i+1));
% 
%      end
%      meanxa_EnKF_CHOLESKY(:,i+1)=mean(Xa_C(:,:,i+1),2);
%      
%   
% %    imagesc(inB)
% %    colorbar
% %    title(sprintf('%i',i))
% %    pause(0.1)
% 
% % plot(i,norm(abs(sum(meanxa_EnKF_CHOLESKY(:,i)-Xreal(:,i)))),'*r');hold on                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
% % plot(i,norm(abs(sum(meanxa_EnKF(:,i)-Xreal(:,i)))),'*k');hold on 
% 
% % plot(i,mean(abs(meanxa_EnKF_CHOLESKY(1,i)-Xreal(1,i)).^2)));hold on 
% end






%% plot the simulation ensemble and the mean
window=10;  %movil mean for the errors  

for i=1:N
plot(squeeze(Xb(1,i,:)))
    hold on
end
plot(meanxb(1,:),'k','LineWidth',3)
figure

% imagesc(B);
% colormap(flipud(hot));colorbar;title('r = 1');xlabel('state number');ylabel('state number')
% caxis([0 0.5]);
% figure
% imagesc(P0);
% colormap(flipud(hot));colorbar;title('r = 1');xlabel('state number');ylabel('state number')
% caxis([0 0.5]);


error_EnKF=norm(abs(sum(meanxa_EnKF(:,:)-Xreal(:,:))))
figure
imagesc(Xreal),title('Truth State')
figure
%imagesc(meanxa_EnKF),title('Analysis State EnKF')
figure
subplot(2,4,1)
plot(Xreal(10,:),'r','LineWidth',2)
hold on
plot(meanxa_EnKF(10,:),':b','LineWidth',2)
hold on
plot(meanxa_EnKF_CHOLESKY(10,:),'--g','LineWidth',2)
%legend({'X truth','Xa ENKF','Xa ENKF MC'})
legend({'X truth','Xa ENKF'})
title('State 10')


subplot(2,4,5)
for i=1:Tsim
Error_ENKF_CHOLESKY(i)=norm(abs(sum(meanxa_EnKF_CHOLESKY(10,i)-Xreal(10,i))));hold on
Error_ENKF(i)=norm(abs(sum(meanxa_EnKF(10,i)-Xreal(10,i))));hold on
set(gca, 'YScale', 'log')
end
plot(movmean(Error_ENKF)
plot(Error_ENKF_CHOLESKY)    
legend({'Xa ENKF','Xa ENKF MC'})

subplot(2,4,2)
plot(Xreal(20,:),'r','LineWidth',2)
hold on
plot(meanxa_EnKF(20,:),':b','LineWidth',2)
hold on
plot(meanxa_EnKF_CHOLESKY(20,:),'--g','LineWidth',2)
legend({'X truth','Xa ENKF','Xa ENKF MC'})
title('State 20')
subplot(2,4,6)
for i=1:Tsim
Error_ENKF_CHOLESKY(i)=norm(abs(sum(meanxa_EnKF_CHOLESKY(20,i)-Xreal(20,i))));hold on
Error_ENKF(i)=norm(abs(sum(meanxa_EnKF(20,i)-Xreal(20,i))));hold on
set(gca, 'YScale', 'log')
end
plot(Error_ENKF)
plot(Error_ENKF_CHOLESKY)
legend({'Xa ENKF','Xa ENKF MC'})

subplot(2,4,3)
plot(Xreal(30,:),'r','LineWidth',2)
hold on
plot(meanxa_EnKF(30,:),':b','LineWidth',2)
hold on
plot(meanxa_EnKF_CHOLESKY(30,:),'--g','LineWidth',2)
legend({'X truth','Xa ENKF','Xa ENKF MC'})
title('State 30')
subplot(2,4,7)
for i=1:Tsim
Error_ENKF_CHOLESKY(i)=norm(abs(sum(meanxa_EnKF_CHOLESKY(30,i)-Xreal(30,i))));hold on
Error_ENKF(i)=norm(abs(sum(meanxa_EnKF(30,i)-Xreal(30,i))));hold on
set(gca, 'YScale', 'log')
end
plot(Error_ENKF)
plot(Error_ENKF_CHOLESKY)
legend({'Xa ENKF','Xa ENKF MC'})

subplot(2,4,4)
plot(Xreal(40,:),'r','LineWidth',2)
hold on
plot(meanxa_EnKF(40,:),':b','LineWidth',2)
hold on
plot(meanxa_EnKF_CHOLESKY(40,:),'--g','LineWidth',2)
legend({'X truth','Xa ENKF','Xa ENKF MC'})
title('State 40')

subplot(2,4,8)
for i=1:Tsim
Error_ENKF_CHOLESKY(i)=norm(abs(sum(meanxa_EnKF_CHOLESKY(40,i)-Xreal(40,i))));hold on
Error_ENKF(i)=norm(abs(sum(meanxa_EnKF(40,i)-Xreal(40,i))));hold on
set(gca, 'YScale', 'log')
end
plot(Error_ENKF)
plot(Error_ENKF_CHOLESKY)
legend({'Xa ENKF','Xa ENKF MC'})


