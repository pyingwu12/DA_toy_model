function pc_rmse(Xtr,Xfe,oer,st,ed,titin,fil)
siz=size(Xfe);
en=siz(2);
%-------------estimate rmse-------------
se=zeros(siz);
for n=1:en
se(:,n,:)=(squeeze(Xfe(:,n,:))-Xtr).^2;
end
% ----
rmse=zeros(siz(1),siz(3));
for i=1:siz(3)
rmse(1,i)=(mean(se(1,:,i))).^.5;
rmse(2,i)=(mean(se(2,:,i))).^.5;
rmse(3,i)=(mean(se(3,:,i))).^.5;
end
% ---mean rmse
merr(1)=mean(rmse(1,:)); merr(2)=mean(rmse(2,:)); merr(3)=mean(rmse(3,:));
%-------------------------------------------
yli=max(max(rmse));
t=(st:ed)*0.01;
% ---
figure('position',[100,100,780,580])
% ---plot x
axes('Units','pixels','position',[80 70 620 450])
  plot(t,rmse(1,:),'r','LineWidth',0.8); hold on
  line([t(1) t(ed-st+1)],[oer(1) oer(1)],'color','k','LineWidth',1.5,'LineStyle','-.'); hold on
  line([t(1) t(ed-st+1)],[merr(1) merr(1)],'color','r','LineWidth',1.5);  
% ---legend
  merrs=num2str(merr(1));
  h=legend('rmse','obs error',['mean rmse = ',merrs]);   
  legend('boxoff'); set(h,'fontsize',12)
% ---  
  set(gca,'fontsize',8,'ylim',[0 yli],'xlim',[t(1) t(ed-st+1)],'LineWidth',0.8)
  xlabel('Analysis time','fontsize',12)
  ylabel('x','fontsize',14,'rot',0)
  tit=['RMSE - X ( ',titin,' )'];    title(tit,'fontsize',12)
  filnam=['rmse_x_',fil];          print('-dpng',filnam,'-r300')
% ---plot y
figure('position',[100,100,780,580])
axes('Units','pixels','position',[80 70 620 450])
  plot(t,rmse(2,:),'r','LineWidth',0.8); hold on
  line([t(1) t(ed-st+1)],[oer(2) oer(2)],'color','k','LineWidth',1.5,'LineStyle','-.'); hold on
  line([t(1) t(ed-st+1)],[merr(2) merr(2)],'color','r','LineWidth',1.5);
% ---legend  
  merrs=num2str(merr(2));
  h=legend('rmse','obs error',['mean rmse = ',merrs]);
  legend('boxoff'); set(h,'fontsize',12)
% --- 
  set(gca,'fontsize',8,'ylim',[0 yli],'xlim',[t(1) t(ed-st+1)],'LineWidth',0.8)
  xlabel('Analysis time','fontsize',12)
  ylabel('y','fontsize',14,'rot',0)
  tit=['RMSE - Y ( ',titin,' )'];   title(tit,'fontsize',12,'rot',0)
  filnam=['rmse_y_',fil];         print('-dpng',filnam,'-r300')
% ---plot z
figure('position',[100,100,780,580])
axes('Units','pixels','position',[80 70 620 450])
  plot(t,rmse(3,:),'r','LineWidth',0.8); 
  line([t(1) t(ed-st+1)],[oer(3) oer(3)],'color','k','LineWidth',1.5,'LineStyle','-.'); hold on
  line([t(1) t(ed-st+1)],[merr(3) merr(3)],'color','r','LineWidth',1.5);
% ---legend  
  merrs=num2str(merr(3));
  h=legend('rmse','obs error',['mean rmse = ',merrs]);
  legend('boxoff'); set(h,'fontsize',12)
% ---
  set(gca,'fontsize',8,'ylim',[0 yli],'xlim',[t(1) t(ed-st+1)],'LineWidth',0.8)
  xlabel('Analysis time','fontsize',12)
  ylabel('z','fontsize',14,'rot',0)
  tit=['RMSE - Z ( ',titin,' )'];    title(tit,'fontsize',12,'rot',0)
  filnam=['rmse_z_',fil];          print('-dpng',filnam,'-r300')
%}
end