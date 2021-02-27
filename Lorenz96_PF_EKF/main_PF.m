%------------------------
% Particle Filter for Lorenz 96 model
% by PY Wu, 2019.10.9
% update(191009): make weighting calculate for all variables at once 
%------------------------
clear; close all
tic
%---model---
global F
F=8;                     % force of model
nvar=40;                 % number of model variables
xint=ones(nvar,1)*F;     % homogenes initial condition
xint(20)=xint(20)+0.1;   % initial perterbation
%---exprimental settings---
dt=0.01;                 % time step
spin=0.3 *365*0.2/dt;    % length of spin up steps (decided by the first number, unit:year) 
anal=0.3 *365*0.2/dt;    % length of DA steps (decided by the first number, unit:year) 
nstep=spin+anal;    % total steps for integration
obstep=0.25 *0.2/dt;     % interval of observation (decided by the first number, unit:day) 
%---
ensize=1000;    % ensemble size
obserr=2^0.5;   % obs error (the same for producing obs and assimilation)
moderr=1;       % initial error of ensemble mean
enpert=4;       % initial ensemble spread
infl=2;         % inflation for SIR
%---Produce truth------------------------------------
Xt=zeros(nvar,nstep);
Xt(:,1)=xint;
for ti=1:nstep-1
  Xt(:,ti+1)=stepit(@lorenz96,Xt(:,ti),dt);    
end
%---Produce observation (at every location)----------------------
obs=zeros(nvar,floor(anal/obstep));
p=0;
for ti=spin+obstep : obstep : spin+anal
   p=p+1; 
   obs(:,p)=Xt(:,ti)+normrnd(0,obserr,[nvar,1]);   
end
%---Perturb initial ensemble-----------------------------------
Xf=zeros(nvar,ensize);  
for vi=1:nvar    
    Xf(vi,:)=Xt(vi,1)+moderr+normrnd(0,enpert,[1,ensize]);    
end
%---ensemble spin up---
for ti=1:spin
   Xf=stepit(@lorenz96en,Xf,dt);  
end

%---noDA ensemble forecast---
Xno=zeros(nvar,ensize,anal);
Xno(:,:,1)=Xf;
for ti=1:anal-1
   Xno(:,:,ti+1)=stepit(@lorenz96en,Xno(:,:,ti),dt);  
end
Xnom=squeeze(mean(Xno,2));
rmsenoda=rmse96(Xnom(1:anal),Xt(:,spin+1:spin+anal));  
rmseobs=rmse96(obs,Xt(:,spin+1:obstep:spin+anal));
% figure; plot(rmsenoda); hold on; plot(1:obstep:anal,)
% set(gca,'Xlim',[1 anal]); title('RMSE of noDA En mean'); legend('noDA','obs error')

%-----PF--------------------------
Xa=zeros(nvar,ensize);  wei=ones(ensize,size(obs,2))/ensize; Neff=zeros(size(obs,2),1); 
Xda=zeros(nvar,ensize,anal); Xdam=zeros(nvar,anal+1);
Xda(:,:,1)=Xf;
p=0;  
for ti=1:anal 
    if mod(ti,obstep)==0
       p=p+1; 
       inno=squeeze(Xda(:,:,ti))-repmat(obs(:,p),1,ensize);
       xw=1/(2*pi*obserr^2)^(0.5) * exp( -mean(inno.^2,1) / (2*obserr^2) ); % use Gaussian distribution for weighting
       if p==1;  wei(:,p)=xw';  else;  wei(:,p)=xw'.*wei(:,p-1); end         
       wei(:,p)=wei(:,p)./sum(wei(:,p));
       Neff(p)=1/sum(wei(:,p).^2);
       %---resampling---
       if Neff(p) < ensize*0.5 
           Xa=resmapling(wei(:,p),Xda(:,:,ti),infl); % the last number is the purterbation of members when they are resampled to the same value
           wei(:,p)=1/ensize;
       else
           Xa=squeeze(Xda(:,:,ti));
       end             
       Xda(:,:,ti+1)=stepit(@lorenz96en,Xa,dt);
       Xdam(:,ti)= sum( Xa.* repmat(wei(:,p)',nvar,1), 2 );  % ensemble weighted mean
    else  % if mod(ti,obstep)==0
       Xda(:,:,ti+1)=stepit(@lorenz96en,Xda(:,:,ti),dt); 
       if p==0;  Xdam(:,ti)=squeeze(mean(Xda(:,:,ti),2));  else;  Xdam(:,ti)= sum( Xda(:,:,ti).* repmat(wei(:,p)',nvar,1), 2 );  end
    end
end
% endPF=toc
%%
%----plot time series of RMSE----
rmse=rmse96(Xdam(:,1:anal),Xt(:,spin+1:spin+anal));  
figure('position',[50 100 1100 400]); 
plot(rmse,'linewidth',1.5); hold on; 
plot(1:obstep:anal,rmseobs,'linewidth',1.5);  
plot(rmsenoda,'k','linewidth',1);
set(gca,'Xlim',[1 anal],'fontsize',14);  xlabel('Step','fontsize',16)
title(['RMSE of ensemble mean (',num2str(ensize),' mem)'],'fontsize',16)
legend('PF','obs','noDA','location','northeastoutside','fontsize',14)
%saveas(gcf,['PF_rmse_m',num2str(ensize),'.png'])
%{
%----plot Hovm?ller diagram of error (with PF)-------------------
figure('position',[100 50 900 500])
errorX=squeeze(  mean(Xda(:,:,1:anal),2)  )-Xt(:,spin+1:spin+anal);
contourf(errorX,15,'linestyle','none')
colorbar
caxis([-10 10])
title('Error of ensemble mean w/ PF','fontsize',14)
xlabel('step','fontsize',14)
%----plot Hovm?ller diagram of error (without PF)-------------------
figure('position',[150 50 900 500])
contourf(squeeze( mean(Xno(:,:,1:anal),2)  )-Xt(:,spin+1:spin+anal),10,'linestyle','none')
colorbar
caxis([-10 10])
title('Error of ensemble mean w/o PF','fontsize',14)
xlabel('step','fontsize',14)
%}
