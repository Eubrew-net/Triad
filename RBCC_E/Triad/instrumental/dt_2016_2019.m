dt=cell(3,1);rs=cell(3,1);
brewer=[157,183,185]
for i=1:3
    dt{i}=[];
    rs{i}=[];
for ano=2016:2019
 
   [dt_,rs_]=read_dt_obs(brewer(i),['/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/',num2str(ano)]);
   dt{i}=[dt{i};dt_];
   rs{i}=[rs{i};rs_
       ];
end

figure(i)
dt_avg=dt{i};
errorbar(dt_avg(:,1),dt_avg(:,32),dt_avg(:,33),'b.')
hold on;
box on;
errorbar(dt_avg(:,1),dt_avg(:,30),dt_avg(:,31),'r.')
%h=ploty(smoothdata(dt_avg(:,[1,30,32]),'rloess',90));
dt_avg_m=meanmonth(dt_avg(:,[1,30,32]));
h=ploty(dt_avg_m.media(:,[1,5,6]))
set(h,'Linewidth',4)
axis 'tight'
datetick('x',12,'keepticks','keeplimits')
ylim(round(median(dt_avg(:,32)))+[-20,20])
title(num2str(brewer(i)))
grid on;
legend('dt low','dt high')
end