function rmse=rmse96(Xfm,Xt)
% Dimension of Xfm (ensemble mean): (vari,time)  
% Dimension of Xt (Truth): (vari,time)
nvar=size(Xt,1);
err= Xfm - Xt;
rmse=(sum(err.^2,1)/nvar).^0.5;
