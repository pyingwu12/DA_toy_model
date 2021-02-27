function Xa=etkf(Xf,Y,en,H,obserr)
% --- Ensemble Transform Kalamn Filter ---
%   input:
%     Xf : model state (background), size=n*k (varibale number*ensemble size)
%     Y : observation, size=p*1
%     en : ensemble size
%     H : observation operator
%     obserr : observation error for each observation varible
%   ouput:
%     Xa : model state (analysis), size=Xf  
%--------------------------------------

% ---set about X
Xf_bar(:,1)=[mean(Xf(1,:)); mean(Xf(2,:)); mean(Xf(3,:))];
Xf_p=zeros(length(Xf(:,1)),en);
for i=1:en
Xf_p(:,i)=(Xf(:,i)-Xf_bar)./(en-1)^0.5;
end
Pf=(Xf_p*Xf_p');
% ---set about obs
Yf=H*Xf; % model value at obs state
Yf_bar(:,1)=[mean(Yf(1,:)); mean(Yf(2,:)); mean(Yf(3,:))];
Yf_p=H*Xf_p;
R=diag(diag(obserr*obserr'));
% ---pert (get Pa)
I=eye(en);
T=sqrtm(inv(I+Yf_p'/R*Yf_p));
Xa_p=Xf_p*T; %***get Pa
% ---mean (get Xa_bar)
K=Pf*H'/(H*Pf*H'+R);
Xa_bar=Xf_bar+K*(Y-Yf_bar);
% ---mean + pert
Xa=zeros(size(Xf));
for i=1:en
Xa(:,i)=Xa_p(:,i).*(en-1)^0.5+Xa_bar;
end

end