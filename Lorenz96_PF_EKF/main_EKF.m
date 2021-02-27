%------------------------
% Extended Kalman filter for Lorenz 96 model
% PY Wu, 2020/09/18
%------------------------
clear;  close all
%---model---
global F
F=8;                    % force of model
nvar=40;                % number of variables
xint=ones(nvar,1)*F;    % initial number
xint(20)=xint(20)+0.1;  % initial perterbation of the model
%---expri---
dt=0.01;                % time step
spin=0.3*365*0.2/dt;     % number of spin up steps
anal=0.5*365*0.2/dt;     % number of DA steps
nstep=spin+anal;   % total steps of integration
obstep=0.25*0.2/dt;     % interval of observation
Nobs=40;        % numbers of obs at one time
obserr=1.141;   % obs error
infl=1.2;        % inflation
%---obs operator, Hx=y---------------
obsint=floor(nvar/Nobs);  p_loc=1:obsint:nvar;
%p_loc=[1 6 7 8 16 18 25 29 38 40]; Nobs = length(p_loc); % for unregular obs
H=zeros(Nobs,nvar);   
for p=1:Nobs
  H(p,p_loc(p))=1;
end
R=(eye(Nobs)*obserr).^2;   % obs error matrix
%---Produce truth and obs------------------------------------
Xt=zeros(nvar,nstep);  Xt(:,1)=xint;
obs=zeros(Nobs,floor(anal/obstep)); p=0;
for ti=1:nstep-1
  Xt(:,ti+1)=stepit(@lorenz96,Xt(:,ti),dt); 
  if ti>spin && ti<= spin+anal && mod(ti-spin,obstep)==0
    p=p+1;
    obs(:,p)=H*Xt(:,ti) + randn(Nobs,1)*obserr;
  end
end
%
%---noDA forecast---
Xno=zeros(nvar,anal);
Xno(:,1)=Xt(:,spin+1)+randn(nvar,1)*0.01;
for ti=1:anal-1
  Xno(:,ti+1)=stepit(@lorenz96,Xno(:,ti),dt);  
end
%%
%-------------------------
Pa=eye(nvar);
Xda=zeros(nvar,anal);  Xda(:,1)=Xno(:,1);
M=zeros(nvar,nvar); %TLM   
delta=0.05;
x=repmat(Xda(:,1),1,nvar);
x_de=x+eye(nvar)*delta;
p=0;  
for ti=1:anal-1 
  if mod(ti,obstep)==0   
    p=p+1;     
    M=(x_de-x)/delta;            % M = ( Mx - M(x+dej) ) / d
    Pf=M*Pa*M';                  % Pf = M Pa M^T     
    d=obs(:,p)-H*Xda(:,ti);      % d = y - Hxb
    K=Pf*H'*inv((R+H*Pf*H'));    % K = Pf H ( R + H Pf H^T )^-1
    Xda(:,ti)=Xda(:,ti)+K*d;     % Xa = Xf + K ( y - Hxb )
    Pa=(eye(nvar)-K*H)*Pf ;      % Pa = ( I - KH ) Pf
    Pa=Pa*infl;                  % inflation        
    x=repmat(Xda(:,ti),1,nvar);
    x_de=x+eye(nvar)*delta;
  end    
  %trPa(ti)=trace(Pa)/nvar;       %  Sum of diagonal elements of Pa
  Xda(:,ti+1)=stepit(@lorenz96,Xda(:,ti),dt);   
  x=stepit(@lorenz96en,x,dt);    % integration for calculation M
  x_de=stepit(@lorenz96en,x_de,dt); 
end
%
%%
plotvari=20; plotlen=1000;
figure('position',[10 50 1000 500]); 
plot(Xt(plotvari,spin:spin+anal),'linewidth',2); hold on
plot(Xno(plotvari,1:anal),'linewidth',1.5)
plot(Xda(plotvari,1:anal),'linewidth',1.5)
plot(obstep:obstep:anal,obs(plotvari,:),'o','linewidth',1.2); %!!! only work when Nobs=40 !!!
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
