function da_plot(t,Xtr,Xf,Xf_s,obs,st,et,sz,titin,fil,ind)
tL=length(t);
figure('position',[50,100,780,580])
% ---------- plot x ----------
axes('Units','pixels','position',[45,385,650,135])
  plot(t(st-sz:tL),Xtr(1,st-sz:tL),'k','LineWidth',1.2); hold on                     % truth
  plot(t(et:tL),Xf(1,et:tL),'color',[0.9 0.2 0.2],'LineWidth',0.8); % da fcst
  plot(t(et:tL),Xf_s(1,:),'color',[0.2 0.2 0.8],'LineWidth',0.8);   % noda fcst
  if ind~=2 && ind~=1 && ind~=3
  plot(t(st:sz:et),obs(1,st:sz:et),'o','Markersize',3,'MarkerFaceColor',[0.8 0.8 1]);  % obs
  h=legend('truth','da','noda','obs');  set(h,'position',[690 500 68 40]);
  end
  plot(t(st-sz:et),Xf(1,st-sz:et),'r','LineWidth',0.8);                     % analysis
% ---setting fig.
  set(gca,'fontsize',8,'ylim',[-20 20])
  ylabel('x','fontsize',12,'rot',0)
  tit=['Trajectory ( ',titin,' )'];   title(tit,'fontsize',12)
% ---------------------------  
% ----------plot y ---------- 
axes('Units','pixels','position',[45,220,650,135])
  plot(t(st-sz:tL),Xtr(2,st-sz:tL),'k','LineWidth',1.2); hold on                     % truth
  plot(t(et:tL),Xf(2,et:tL),'color',[0.9 0.2 0.2],'LineWidth',0.8);
  plot(t(et:tL),Xf_s(2,:),'color',[0.2 0.2 0.8],'LineWidth',0.8);  
  if ind~=4 && ind~=1 && ind~=5
  plot(t(st:sz:et),obs(2,st:sz:et),'o','Markersize',3,'MarkerFaceColor',[0.8 0.8 1]);  %obs 
  h=legend('truth','da','noda','obs');   set(h,'position',[690 500 68 40]);
  end
  plot(t(st-sz:et),Xf(2,st-sz:et),'r','LineWidth',0.8);                     % analysis
% ---setting fig.
  set(gca,'fontsize',8,'ylim',[-30 30])
  ylabel('y','fontsize',12,'rot',0)
% --------------------------- 
% -----------plot z ---------
axes('Units','pixels','position',[45,50,650,135])
  plot(t(st-sz:tL),Xtr(3,st-sz:tL),'k','LineWidth',1.2); hold on                     % truth
  plot(t(et:tL),Xf(3,et:tL),'color',[0.9 0.2 0.2],'LineWidth',0.8);
  plot(t(et:tL),Xf_s(3,:),'color',[0.2 0.2 0.8],'LineWidth',0.8);  
  if ind~=4 && ind~=2 && ind~=6
  plot(t(st:sz:et),obs(3,st:sz:et),'o','Markersize',3,'MarkerFaceColor',[0.8 0.8 1]);  %obs
  h=legend('truth','da','noda','obs');   set(h,'position',[690 500 68 40]);
  end
  plot(t(st-sz:et),Xf(3,st-sz:et),'r','LineWidth',0.8);                     % analysis
% ---setting fig.  
  set(gca,'fontsize',8,'ylim',[0 50])
  xlabel('time','fontsize',12)
  ylabel('z','fontsize',12,'rot',0)
%----------------------------
  filnam=['trajectory_',fil];   print('-dpng',filnam,'-r600')
%}
end