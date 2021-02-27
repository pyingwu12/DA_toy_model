function [t,Y]=odefun(dydt,t0,y0,dt,n,method)
%  hw10_ode_fun(dydt,t0,y0,dt,n,method)
%  dydt ode eqution (function)
%  t0 initial t (scalar)
%  y0 initial value (vector)
%  dt step size
%  n number for time step
%  method (1)euler (2)leapfrog (3)rk-4

t=t0+(0:n)*dt;
t=t';
Y=zeros(length(y0),n+1);
Y(:,1)=y0;
   switch(method)
   case(1)      
      for i=1:n
          f1=feval(dydt,t(i),Y(:,i));
          Y(:,i+1)=Y(:,i)+dt*f1;
      end
   case(2)
      Y(:,2)=Y(:,1)+dydt(t(1),Y(:,1))*dt;
      for i=2:n
        Y(:,i+1)=Y(:,i-1)+dydt(t(i),Y(:,i))*dt*2;
      end        
   case(3)
       for i=1:n
         k1=dt*dydt(t(i),Y(:,i));
         k2=dt*dydt(t(i)+dt/2, Y(:,i)+k1/2);
         k3=dt*dydt(t(i)+dt/2, Y(:,i)+k2/2);
         k4=dt*dydt(t(i)+dt,Y(:,i)+k3);
         Y(:,i+1)=Y(:,i)+(k1+2*k2+2*k3+k4)/6;
       end
   end
t=t';
end


