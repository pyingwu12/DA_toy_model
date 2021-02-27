function xout = lorenz96en(xin)
% Lorenz 96 model, calculate several members at the same time
% input: <xin> should be a matrix with size Nvar*ensize, 
% which Nvar=number of variables, ensize=number of members
% PY Wu, 2020/09/18

global F
[N, ensize]=size(xin);
xout=zeros(N,ensize);

%first the 3 edge cases: i=1,2,N
xout(1,:)=xin(N,:).*(xin(2,:)-xin(N-1,:))-xin(1,:);
xout(2,:)=xin(1,:).*(xin(3,:)-xin(N,:))-xin(2,:);
xout(N,:)=xin(N-1,:).*(xin(1,:)-xin(N-2,:))-xin(N,:);

%then the general case
for i=3:N-1
 xout(i,:)=xin(i-1,:).*(xin(i+1,:)-xin(i-2,:))-xin(i,:);
end

xout=xout+F;