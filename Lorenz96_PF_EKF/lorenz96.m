function xout = lorenz96(xin)
% xin should be a vector

global F
N=size(xin,1);
xout=zeros(N,1);

%first the 3 edge cases: i=1,2,N
xout(1)=xin(N)*(xin(2)-xin(N-1))-xin(1);
xout(2)=xin(1)*(xin(3)-xin(N))-xin(2);
xout(N)=xin(N-1)*(xin(1)-xin(N-2))-xin(N);

%then the general case
for i=3:N-1
 xout(i)=xin(i-1)*(xin(i+1)-xin(i-2))-xin(i);
end

xout=xout+F;
