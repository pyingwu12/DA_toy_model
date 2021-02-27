function au0=forward_adj(au1,u1,dt,dx,gamma,xnu)

Ndim=length(au1);
au0=zeros(Ndim,1);
aA=0; aD=0; aS=0; 

for i=1:Ndim %space    
  ip=i+1;  if i==Ndim; ip=1; end
  ipp=i+2; if i>=Ndim-1; ipp=i+2-Ndim; end
  im=i-1; if i==1; im=Ndim; end
  imm=i-2; if i<=2; imm=i-2+Ndim; end
  %---
  
  au0(i)=au0(i)+au1(i);
  aA=aA + dt*au1(i);
  aD=aD + dt*au1(i);
  aS=aS + dt*au1(i);
  au1(i)=0;
  
  au0(ip)= au0(ip) + xnu/dx^2 * aS;
  au0(i) = au0(i) - 2*xnu/dx^2 * aS;
  au0(im)= au0(im) + xnu/dx^2 * aS; 
  aS=0;
  
  au0(ipp) = au0(ipp) - gamma^2/(2*dx^3) * aD ;
  au0(ip) = au0(ip) + 2*gamma^2/(2*dx^3) * aD ;
  au0(im) = au0(im) - 2*gamma^2/(2*dx^3) * aD ;
  au0(imm) = au0(imm) + gamma^2/(2*dx^3) * aD ;
  aD=0;
  
  au0(ip)= au0(ip) + (-2*u1(ip)-u1(i))* aA / (6*dx);
  au0(i)= au0(i) - (u1(ip)-u1(im))* aA / (6*dx);
  au0(im)= au0(im) + (u1(i)+2*u1(im))* aA / (6*dx);
  aA=0;
  
end    %Ndim  