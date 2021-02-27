function u1=forward(u0,dt,dx,gamma,xnu)
 
Ndim=length(u0);
u1=zeros(Ndim,1);
for i=1:Ndim %space    
  ip=i+1;  if i==Ndim; ip=1; end
  ipp=i+2; if i>=Ndim-1; ipp=i+2-Ndim; end
  im=i-1; if i==1; im=Ndim; end
  imm=i-2; if i<=2; imm=i-2+Ndim; end
  %---
  A = -(u0(ip)+u0(i)+u0(im))*(u0(ip)-u0(im)) / (6*dx) ;
  D = -gamma^2 * (u0(ipp)-2*u0(ip)+2*u0(im)-u0(imm)) / (2*dx^3);
  S = xnu * (u0(ip)-2*u0(i)+u0(im)) / dx^2;     

  u1(i) = u0(i) + dt * (A + D + S);  
    
end %Ndim  