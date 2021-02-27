function Yout=stepit(dydt,yin,dt)
  k1=dt*dydt(yin);
  k2=dt*dydt(yin+k1/2);
  k3=dt*dydt(yin+k2/2);
  k4=dt*dydt(yin+k3);
  Yout=yin+(k1+2*k2+2*k3+k4)/6;
end