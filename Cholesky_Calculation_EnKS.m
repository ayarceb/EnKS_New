%  Code to calculate the snapshots and Cholesky decomposition matrix of the
%  snapshots form the model Lorenz96
%  Andrés Yarce Botero 6/03/2020

clear all
close all
clc

%%  Input paramater configurations

Tsim=500;        % Simulation time
dt=0.001;         % Step length
m=30;            % Number of observations
n=40;            % State number
F=8;             % Forcing factor model
sigma=1e-2;
R=sigma^2*eye(m);       % Observation Covariance matrix 
H = eye(n,n);           % Observation operator matrix
H = H(randperm(n,m),:);
imagesc(H)    %plot measurenment states to know which states are being observed
xlabel('State number');
ylabel('Observation');

%===Generate real state===
x0=5*rand(n,1);
[Xreal]=Lorenz_96(Tsim,dt,x0,F);            % Create vector true


%==Number of Ensembles==
N=100;                    % Number of ensembles
Xb=zeros(n,N,Tsim);       % Background ensemble initialization
Xb(:,:,1)=1*rand(n,N);    
Xa=Xb;
meanXB=zeros(n,Tsim);      % Matrix to put the mean of the background ensemble


for i=1:Tsim-1
    i
    %===== Forecast Step=====
    for en=1:N
        [Xb(:,en,i+1)]=Lorenz_96_one_step(1,dt,squeeze(Xa(:,en,i)),F);
     end 
%     meanxb=mean(Xb(:,:,i),2);
% %     meanXB(:,i)=meanxb;
%     L=Xb(:,:,i)-meanxb;     % Deviation matrix
%     P0=((1/N-1)*(L*L'));    % Background Covariance
%     B=P0; 
%    % ===== Analysis Step=====
%     K=B*H'*pinv(H*B*H'+R);
%    for en=1:N        
%        C(:,i)= Y(:,i)+sigma*normrnd(0,sqrt(1),[1,m])';  % Atempt testing other normal noise generator
%     
%        %C(:,i)= Y(:,i+1)+sigma*randn(m,0.1);
%         Xa(:,en,i+1)=Xb(:,en,i+1)+K*(C(:,i)-H*Xb(:,en,i+1));
%     end
%     meanxa_EnKF(:,i+1)=mean(Xa(:,:,i+1),2);
end