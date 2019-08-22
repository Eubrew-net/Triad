dt=cell(3,1);rs=cell(3,1);
brewer=[157,183,185]
read_config_
for i=1:3
    dt{i}=[];
    rs{i}=[];
% for ano=2016:2019
%  
%    [dt_,rs_]=read_dt_obs(brewer(i),['/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/',num2str(ano)]);
%    dt{i}=[dt{i};dt_];
%    rs{i}=[rs{i};rs_
%        ];
% end

load('dt2016_2019')
load('rs2016_2019')

 dt_ev=meanperiods(dt{i}, events{i}); 
 rs_ev=meanperiods(rs{i}, events{i}); 

  data=[events{i}.dates round(dt_ev.m(:,30),1) round(dt_ev.std(:,30),2)  1e9*op_cfg{i}{9,:}' ...
      round(dt_ev.m(:,32),1)     round(dt_ev.std(:,32),2)  1e9*alt_cfg{i}{9,:}' (dt_ev.N(:,end))];
    dt_evt{i}=array2table(data,'VariableNames',{'Date','DT_h','DT_h_std','DT_op_ref','DT_l','DT_l_std','DT_alt_ref','N'},'RowNames',varname(str2name(strrep(events{i}.labels,'"',''))));
    dt_evt{i}.Fecha=datetime(datestr(dt_evt{i}.Date));
    dt_evt{i}=timetable2table(table2timetable(dt_evt{i}));
    writetable(dt_evt{i},'IzoTriad_2016_2019.xls','Sheet',strrep('dt_157','157',num2str(brewer(i))),'WriteRowNames',true)
 
figure(i)
dt_avg=dt{i};
errorbar(dt_avg(:,1),dt_avg(:,32),dt_avg(:,33),'b.')
hold on;
box on;
errorbar(dt_avg(:,1),dt_avg(:,30),dt_avg(:,31),'r.')
%h=ploty(smoothdata(dt_avg(:,[1,30,32]),'rloess',90));
dt_avg_m=meanmonth(dt_avg(:,[1,30,32]));
h=ploty(dt_avg_m.media(:,[1,5,6]));
set(h,'Linewidth',4);
set(h(2),'Color', [0 1 0]);
set(h(1),'Color', [0 0 0]);


axis 'tight'
datetick('x',12,'keepticks','keeplimits')
ylim(round(median(dt_avg(:,32)))+[-20,20])
title(num2str(brewer(i)))
grid on;
legend('dt low','dt high')
end