function [Xtr, Xien, obs, t]=makedata(en,obserr,backerr,n)
%---------------Lorenz & truth state--------------
% set parameter
dt=0.01;           %time step
Yi=[10 0 38];      %initial number
%function
r=28; s=10; b=8/3; %coefficient
dy=@(t,y) [-s*y(1)+s*y(2) ;r*y(1)-y(2)-y(1)*y(3) ;-b*y(3)+y(1)*y(2)];
% --- make truth state
[t,Xtr]=feval('odefun',dy,0,Yi,dt,n,3); 
xt=Xtr(1,:); yt=Xtr(2,:); zt=Xtr(3,:);
%-------------------------------------------------
%-----------make obs & initial condition----------
% --- make obs
  obs(1,:)=xt+(normrnd(0,obserr(1),[1,n+1]));
  obs(2,:)=yt+(normrnd(0,obserr(2),[1,n+1]));
  obs(3,:)=zt+(normrnd(0,obserr(3),[1,n+1]));
% --- make ensemble initial condition
  Xien(1,:)=xt(1)+5+normrnd(0,backerr,[en,1]);
  Xien(2,:)=yt(1)+5+normrnd(0,backerr,[en,1]);  
  Xien(3,:)=zt(1)+5+normrnd(0,backerr,[en,1]);  
%---------------------------------------------------
end