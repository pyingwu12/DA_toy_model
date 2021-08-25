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
rho=1; % inflation, change this to an input variable of the function if the inflation is needed
%---set about X
Xf_bar(:,1)=mean(Xf,2);
Xf_p=zeros(length(Xf(:,1)),en);
for i=1:en
 Xf_p(:,i)=(Xf(:,i)-Xf_bar)./(en-1)^0.5;
end
%---model state at obs space
Yf_bar=H*Xf_bar; 
Yf_p=H*Xf_p;
%--obs error maxtrix
R=diag(diag(obserr*obserr'));
%---update ensemble perturbation (get Pa)
I=eye(en);
Pa_hat=inv(I/rho+Yf_p'/R*Yf_p);
W=sqrtm(Pa_hat);
Xa_p=Xf_p*W; 
%---update ensemble mean (get Xa_bar)
w_bar=Pa_hat*Yf_p'/R*(Y-Yf_bar);
Xa_bar=Xf_bar+Xf_p*w_bar;
% ---mean + pert
Xa=zeros(size(Xf));
for i=1:en
Xa(:,i)=Xa_bar+Xa_p(:,i).*(en-1)^0.5;
end

end