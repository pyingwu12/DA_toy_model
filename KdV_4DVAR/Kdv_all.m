%Kdv
clear
close all
%=============basic setting================
%---model
Ndim = 40;
dx = 10.0;
dt = 0.5;
KTmax=1200; % total integration step
gamma = 20.0;
xnu = 0.2;

%---obs
Nobs = 8;                       % obs numbers in space
s_obsint = floor(Ndim / Nobs);  % obs space interval
t_obsint = 15;                  % obs time interval (step)
robs=0.5;                       % obs error
%---obs operator, Hx=y
Hmat=zeros(Nobs,Ndim); 
p_loc=1:s_obsint:Ndim;
%p_loc=[1 6 7 8 15 16 18 25 27 29 33 38 40]; Nobs = length(p_loc); %irregular obs location
for p=1:Nobs
Hmat(p,p_loc(p))=1;
end
%%
%=============creat true state and obs================
%---initial condition---
x0=200; alpha=16;
ut=zeros(Ndim,KTmax+1);
for i=1:Ndim
ut(i,1)=alpha/(exp((i*dx-x0)/(48/alpha)^0.5/gamma)+exp(-(i*dx-x0)/(48/alpha)^0.5/gamma))^2;
end
%---make true state and obs----
p=0; yobs=zeros(Nobs,KTmax/t_obsint);
for t=1:KTmax     
  if t==1
   ut(:,t+1)=forward(ut(:,t),dt,dx,gamma,xnu);  
  else
   ut(:,t+1)=leapfrog(ut(:,t-1),ut(:,t),dt,dx,gamma,xnu);
  end
  if mod(t,t_obsint)==0
     p=p+1;
     yobs(:,p)=Hmat * ut(:,t+1) + robs * randn(Nobs,1);
  end
end %time
%save(['true_obs_',num2str(Nobs),'_',num2str(t_obsint),'_',num2str(robs),'.mat'],'ut','yobs')
%%
%=============no DA================
%---different initial condiction from true
x0=10; b=110; a=5;
uf=zeros(Ndim,KTmax+1);
for i=1:Ndim
uf(i,1)=a*exp( -((i*dx-x0)/b)^2 );
end
%---free forecast
for t=1:KTmax        
  if t==1
   uf(:,t+1)=forward(uf(:,t),dt,dx,gamma,xnu);  
  else
   uf(:,t+1)=leapfrog(uf(:,t-1),uf(:,t),dt,dx,gamma,xnu);
  end      
end %time
%%
%=============DA================
% load(['true_obs_',num2str(Nobs),'_',num2str(t_obsint),'_',num2str(robs),'.mat']);
%---DA parameters
KTwin = 150;                   % Assimilation window
robs = 0.5;                    % obs error
sigma_f=2;                     % background error
Niter = 8000;                  % Number of iterations
epsilon = 0.0005;             % mutiply to dJ/dx
ecri = 1.0e-6;                 % criteria of iterations

%---background error covariance
L=2;  B=zeros(Ndim,Ndim);
for i=1:Ndim
  for j=1:Ndim
    jj=j-i;
    if jj+Ndim/2>Ndim; jj=jj-Ndim; end
    if jj+Ndim/2<1; jj=jj+Ndim; end
    B(i,j)=sigma_f^2*exp(-abs(jj)./(2*L^2));
  end
end
B=B+eye(Ndim)*0.00001;
%B=sigma_f^2*eye(Ndim);
% figure; imagesc(B); colorbar

%---run DA
uest=zeros(Ndim,KTmax+1);  uest(:,1)=uf(:,1);
disp('---DA START---')
for kt=1:KTwin:KTmax
%for kt=1:KTwin:KTwin*2
  for m=1:Niter
    Jcost = 0; 
 %---forward----------------------------------------     
    for k=1:KTwin     
      ktk=kt+k-1;
      if ktk>KTmax; break; end
      if k==1
       uest(:,ktk+1)=forward(uest(:,ktk),dt,dx,gamma,xnu);  
      else
       uest(:,ktk+1)=leapfrog(uest(:,ktk-1),uest(:,ktk),dt,dx,gamma,xnu);
      end     
      if mod(ktk,t_obsint)==0
        inno=yobs(:,fix(ktk/t_obsint)) - Hmat * uest(:,ktk+1);      
        Jcost=Jcost+(inno'*inno)/robs^2;        
        %Jcost_o(m,1)=(inno'*inno)/robs^2;      
      end    
    end  %k=1:KTwin
    
 %---adjoint---------------------------------------- 
    au2=zeros(Ndim,1); au0T=0; % adjoint variables
    for k=KTwin:-1:1       
      ktk=kt-1+k; 
      if ktk>KTmax; continue; end
      if mod(ktk,t_obsint)==0
        inno2=yobs(:,fix(ktk/t_obsint))- Hmat * uest(:,ktk+1);   
        au2 = au2 +  Hmat' * inno2 ./ robs^2;           
      end  %if mod  
      %-----------
      if k==1
        au2=forward_adj(au2,uest(:,ktk),dt,dx,gamma,xnu) ;  
      else
        [au0, au2]=leapfrog_adj(au2,uest(:,ktk),dt,dx,gamma,xnu);          
      end
      au2=au2+au0T;
      au0T=au0;      
    end  %k=1:KTwin
   
    delx=epsilon * (au2 - inv(B) * uest(:,kt)) ;   
    uest(:,kt)=uest(:,kt)+delx;
    
    if sum(delx.^2)<ecri; break;  end
  end %iterations
  disp([num2str(kt),'th step done, ',num2str(m),' iter']);
  
 %---forcast from analysis---
  for k=1:KTwin         
    ktk=kt-1+k;
    if ktk>KTmax; break; end
    if k==1
     uest(:,ktk+1)=forward(uest(:,ktk),dt,dx,gamma,xnu);  
    else
     uest(:,ktk+1)=leapfrog(uest(:,ktk-1),uest(:,ktk),dt,dx,gamma,xnu);
    end 
  end 
  
end  %1:KTwin:KTmax-1  
disp('---DA DONE---')

%
%============= PLOT ================
exp_setting_tit=['KTwin= ',num2str(KTwin),', robs=',num2str(robs),', Nobs=',num2str(Nobs),', t_obsint=',num2str(t_obsint),...
    ' sigma_f=',num2str(sigma_f),', L=',num2str(L)];
exp_setting_file=['w',num2str(KTwin),'_ro',num2str(robs*10,'%.2d'),'_No',num2str(Nobs),'_to',num2str(t_obsint),'_rf',...
    num2str(sigma_f),'_L',num2str(L)];
filename='4DVARmod';
%%
%---True, forecasr, and analysis results
col_range=[0 6];
hf1=figure('position',[50 50 1100 600]);
subplot(1,3,1)
contourf(ut',20,'linestyle','none')
colorbar; caxis(col_range)
set(gca,'LineWidth',1.2,'FontSize',13)
ylabel('Time step','FontSize',16)
xlabel('X (model grid)')
title('True')

subplot(1,3,2)
contourf(uf',20,'linestyle','none')
colorbar; caxis(col_range)
set(gca,'LineWidth',1.2,'FontSize',13)
xlabel('X (model grid)')
title('Forecast')

subplot(1,3,3)
contourf(uest',20,'linestyle','none')
colorbar; caxis(col_range)
set(gca,'LineWidth',1.2,'FontSize',13)
xlabel('X (model grid)')
title('Analysis')

sgtitle(exp_setting_tit,'Interpreter','none','fontsize',14)
% outfile=[filename,'_',exp_setting_file]; 
% print(hf1,'-dpng',[outfile,'.png'])
%%
%---RMSE time series
rmse_a=(sum((uest-ut).^2,1)/Ndim).^0.5;
rmse_f=(sum((uf-ut).^2,1)/Ndim).^0.5;

hf2=figure('position',[100 50 1000 600]);
plot(1:KTmax+1,rmse_f,1:KTmax+1,rmse_a,'LineWidth',2.5);
legend('forecast','analysis','fontsize',16)
set(gca,'LineWidth',1.2,'FontSize',16,'Xlim',[1 KTmax+1])
xlabel('Time step','FontSize',18)
title({'RMSE',exp_setting_tit},'Interpreter','none','fontsize',15)
% outfile=[filename,'_rmse_',exp_setting_file]; 
% print(hf2,'-dpng',[outfile,'.png'])
%%
%---plot states at initial time of a DA window (NTWIN)
hf3=figure('position',[200 50 850 600]);
NTWIN=0;  t=1+NTWIN*KTwin;
plot(uf(:,t),'LineWidth',2.5); hold on
plot(uest(:,t),'LineWidth',2.5)
plot(ut(:,t),'LineWidth',2.5)
plot(p_loc,yobs(:,fix(t/t_obsint)+1),'O','LineWidth',2)
legend('background','analysis','truth','obs')
set(gca,'LineWidth',1.2,'FontSize',16)
xlabel('X (model grid)')
title({['States at t=',num2str(t)],exp_setting_tit},'Interpreter','none','fontsize',14)
%set(gca,'YLim',[-1 5])
% outfile=[filename,'_t',num2str(t),'_',exp_setting_file]; 
% print(hf3,'-dpng',[outfile,'.png'])
