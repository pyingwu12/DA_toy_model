function x1=stepit(dxdt,x,dt)
  k1=dt*dxdt(x);
  k2=dt*dxdt(x+k1/2);
  k3=dt*dxdt(x+k2/2);
  k4=dt*dxdt(x+k3);
  x1=x+(k1+2*k2+2*k3+k4)/6;
end