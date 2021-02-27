function u2=leapfrog(u0,u1,dt,dx,gamma,xnu)
 
Ndim=length(u0);
u2=zeros(Ndim,1);
for i=1:Ndim %space    
  ip=i+1;  if i==Ndim; ip=1; end
  ipp=i+2; if i>=Ndim-1; ipp=i+2-Ndim; end
  im=i-1; if i==1; im=Ndim; end
  imm=i-2; if i<=2; imm=i-2+Ndim; end
  %---
  A = -(u1(ip)+u1(i)+u1(im))*(u1(ip)-u1(im)) / (6*dx) ;
  D = -gamma^2 * (u1(ipp)-2*u1(ip)+2*u1(im)-u1(imm)) / (2*dx^3);
  S = xnu * (u0(ip)-2*u0(i)+u0(im)) / dx^2;     

  u2(i) = u0(i) + 2 * dt * (A + D + S);  
    
end %Ndim  