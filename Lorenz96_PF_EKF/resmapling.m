function Xa=resmapling(wei_s,Xi,infla)
% Dimension of Xi: (nvar,nsize,ti)

Xi=squeeze(Xi);
ensize=size(Xi,2); nvar=size(Xi,1);
Xa=zeros(nvar,ensize); 
sir=zeros(1,ensize);
% using random number to find important member (SIR)
for j=1:ensize
   rn=rand;      
   sumw=0;
   for mi=1:ensize
      sumw=sumw+wei_s(mi);
      if sumw>=rn;   sir(mi)=sir(mi)+1;   break;  end                         
   end
end
% resample the ensemble by sir found above
respl=find(sir~=0); 
m=0;
for i=1:length(respl)
    for k=1:sir(respl(i))                
      m=m+1;
      Xa(:,m)=Xi(:,respl(i))+normrnd(0,infla);
    end
end
             


