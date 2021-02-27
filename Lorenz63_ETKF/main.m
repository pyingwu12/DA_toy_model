%--------------------------------------------------------------
%     Lorenz 3-variable system with data assimilation
%     DA scheme: ETKF
%     created by Pin-Ying, Wu 2014 (for DA project)
%--------------------------------------------------------------
%
clear
close all
%------------------set parameter------------------
obserr=[1.414;1.414;1.414]; % observation error
backerr=2;                  % background error
anatime=1000;            % analysis time(step)
sptime=500;              % spin up time(step)
%-------------------------------------------------
%------------------data setting-------------------
disp('--------------------STAR---------------------')
% ---new data or not
dataind=input('Use old data(for truth and I.C.)? (1)Y  (2)N : ');
 % ---download old data
if dataind==1; datanam=input('old data file name ? (keyin 1 for defult: data_en10.mat ) : ','s');
  if str2double(datanam)==1; load('data_en10.mat');  else; load([datanam,'.mat']);  end 
  load('obs.mat')
  en=length(Xien(1,:));
 % ---make new data
elseif dataind==2
  en=input('ensemble size = ');
  n=input(['total step (> ',num2str(anatime+sptime),' ) = ']); n=n-1;
  % ---
  [Xtr, Xien, obs, t]=feval('makedata',en,obserr,backerr,n);
  % ---save new data or not
  savind=input('Save new data ? (1)Y (2)N  : ');
  if savind==1; datanam=input('svae file name? : ','s');
     enm=num2str(en); save([datanam,'.mat'],'Xtr','Xien','t');
  end   
  savoind=input('Save new obs? (cover the original one!) (1)Y(use new obs) (2)N(use old obs)  : ');  
  if savoind==1; save('obs.mat','obs'); end 
  load('obs.mat')
end
xt=Xtr(1,:); yt=Xtr(2,:); zt=Xtr(3,:);
tL=length(xt);
%-------------------------------------------------
%--------------experiment setting-----------------
disp(' --- ')
disp('Which obs you want to assimilation?')
daind=input('x:4  y:2  z:1 (ex.y+z: 3, x: 4 ) : '); 
  if daind==4; H=[1 0 0; 0 0 0; 0 0 0];
  elseif daind==2;  H=[0 0 0; 0 1 0; 0 0 0];
  elseif daind==1;  H=[0 0 0; 0 0 0; 0 0 1];
  elseif daind==6;  H=[1 0 0; 0 1 0; 0 0 0];
  elseif daind==5;  H=[1 0 0; 0 0 0; 0 0 1];
  elseif daind==3;  H=[0 0 0; 0 1 0; 0 0 1];
  elseif daind==7;  H=[1 0 0; 0 1 0; 0 0 1]; 
  end
% ---
expind=input('Assimilation step size? (1)defult:8  (2)key in by yourself : ');
if expind==1;      dasz=8; 
elseif expind==2;  dasz=input('assimilation step size = ');
else; disp('--assimilation step size = defult--'); dasz=8;
end
anacyle=floor((anatime)/dasz);
% ---
titin=input('Figure title ? : ','s');
namein=input('Figure file name ? : ','s');
% ---
disp('Function start.....')
%-------------------------------------------------
%---------------------Lorenz-----------------------
% set parameter
r=28; s=10; b=8/3; % coefficient
dt=0.01;           % time step
dy=@(y) [-s*y(1)+s*y(2) ;r*y(1)-y(2)-y(1)*y(3) ;-b*y(3)+y(1)*y(2)];
%--------------------------------------------------
%----------forecast and data assimilation----------
% -------spin up-------
Xfe(:,:,1)=Xien;
ind=0;
for nsp=1:sptime
  ind=ind+1;
  for n=1:en
   Xfe(:,n,ind+1)=feval('stepit',dy,Xfe(:,n,ind),dt);
  end  
end
% -------analysis cycle-------
dastime=ind+dasz+1;
for ana = 1 : anacyle
  for i=1:dasz
  ind=ind+1;
  for n=1:en
   Xfe(:,n,ind+1)=feval('stepit',dy,Xfe(:,n,ind),dt);
  end  
  end
  % --- assmilation
  Xfe(:,:,ind+1)=feval('etkf',Xfe(:,:,ind+1),obs(:,ind+1),en,H,obserr);
end
if ind+1 < anatime+sptime+1
  for i=1:anatime+sptime-ind
  ind=ind+1;
  for n=1:en
  Xfe(:,n,ind+1)=feval('stepit',dy,Xfe(:,n,ind),dt);
  end  
  end     
end
daetime=ind+1;
% -------forecast run-------
%
frind=input('(1)Ensemble forecast or (2)single forecast ? : ');
if frind==1
%--- ensemble fcst
 for i=1:tL-ind-1
  ind=ind+1;
  for n=1:en
  Xfe(:,n,ind+1)=feval('stepit',dy,Xfe(:,n,ind),dt);
  end  
 end 
 enend=tL;   
%
elseif frind==2
% --- single fcst
 Xf=zeros(3,tL);
 Xf(:,daetime)=[mean(Xfe(1,:,daetime)) mean(Xfe(2,:,daetime)) mean(Xfe(3,:,daetime))];
 for i=1:tL-ind-1
  ind=ind+1;
  Xf(:,ind+1)=feval('stepit',dy,Xf(:,ind),dt);
  for n=1:en
  Xfe(:,n,ind+1)=Xf(:,ind+1); % expand Xf to ensemble for estimate rmse  
  end
 end
 enend=daetime-1;
end
% --- set Xf for plot trajectory
 for i=1:enend
  Xf(1,i)=mean(Xfe(1,:,i));   Xf(2,i)=mean(Xfe(2,:,i));   Xf(3,i)=mean(Xfe(3,:,i));
 end 
%--------------------------------------------------
%-------------------noda forecast------------------
Xf_s(:,1)=[obs(1,daetime) obs(2,daetime) obs(3,daetime)];
 for i=1:tL-daetime
    Xf_s(:,i+1)=feval('stepit',dy,Xf_s (:,i),dt);
 end
%--------------------------------------------------
%---------------------plot-------------------------
disp('Plot.....')
da_plot(t,Xtr,Xf,Xf_s,obs,dastime,daetime,dasz,titin,namein,daind)  % truth and forecast
% --- rmse
pc_rmse(Xtr(:,dastime:daetime),Xfe(:,:,dastime:daetime),obserr,dastime,daetime,titin,namein)
%--------------------------------------------------
%}
disp('--------------------END---------------------')