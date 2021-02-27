%------------------------
% Kalman filter for Lorenz 96 model
% by PYWu, 2020.09.09
%------------------------
clear; close all
%---model---
global F
F=8;                     % force of model
nvar=40;                 % number of variables
xint=ones(nvar,1)*F;     % initial number
xint(20)=xint(20)+0.1;   % initial perterbation of the model
%---expri---
dt=0.01;                 % time step
spin=0.3*365*0.2/dt;     % number of spin up steps
anal=0.5*365*0.2/dt;     % number of DA steps
nstep=spin+anal;    % total steps of integration
obstep=0.25*0.2/dt;      % interval of observation
obserr=1.141;   % obs error
P=eye(nvar);    % background error covariance (set to unit diagonal matrix)
R=(eye(nvar)*obserr).^2; % obs error covariance
%---Produce truth------------------------------------
Xt=zeros(nvar,nstep);
Xt(:,1)=xint;
for ti=1:nstep-1
  Xt(:,ti+1)=stepit(@lorenz96,Xt(:,ti),dt);    
end
%---Produce observation------------------------------------
obs=zeros(nvar,floor(anal/obstep));
p=0; H=eye(nvar);
for ti=spin+obstep:obstep:spin+anal
   p=p+1; 
   obs(:,p)=H*Xt(:,ti)+normrnd(0,obserr,[nvar,1]);   
end
%---noDA forecast---
Xno=zeros(nvar,anal);
Xno(:,1)=Xt(:,spin+1)+randn(nvar,1)*0.01;
for ti=1:anal-1
   Xno(:,ti+1)=stepit(@lorenz96,Xno(:,ti),dt);  
end
%--- DA ----------------------------
Xda=zeros(nvar,anal); 
Xda(:,1)=Xno(:,1);
p=0;  
for ti=1:anal-1    
  if mod(ti,obstep)==0
    p=p+1;        
    d=obs(:,p)-H*Xda(:,ti); % innovation
    K=P*H'*inv((R+H*P*H')); % Kalman gain
    Xda(:,ti)=Xda(:,ti)+K*d;
    Xda(:,ti+1)=stepit(@lorenz96,Xda(:,ti),dt); 
  else
    Xda(:,ti+1)=stepit(@lorenz96,Xda(:,ti),dt); 
  end
end
%
%%
plotvari=20; plotlen=1000;
figure('position',[10 50 1000 500]); 
plot(Xt(plotvari,spin+1:spin+anal),'linewidth',2); hold on
plot(Xno(plotvari,1:anal),'linewidth',1.5)
plot(Xda(plotvari,1:anal),'linewidth',1.5)
plot(obstep:obstep:anal,obs(plotvari,:),'o','linewidth',1.2)
legend('Truth','noDA','DA','obs','location','best')
set(gca,'Xlim',[0 plotlen],'fontsize',16,'linewidth',1)
title(['Time series of ',num2str(plotvari),'th variable'],'fontsize',18)
xlabel('Time step')
%%
errno=Xno-Xt(:,spin+1:spin+anal);  rmse_no=(mean(errno.^2,1)).^0.5;  %rmse of noDA
errda=Xda-Xt(:,spin+1:spin+anal);  rmse_da=(mean(errda.^2,1)).^0.5;  %rmse of DA
%---
plotlen=2000;
figure('position',[10 50 1000 500]); 
plot(rmse_no,'linewidth',2); hold on
plot(rmse_da,'linewidth',2)
legend('noDA','DA','location','best')
set(gca,'Xlim',[0 plotlen],'fontsize',16,'linewidth',1)
title('RMSE','fontsize',18)
xlabel('Time step')