close all;
clear all;
ap=cell(3,1);
save('ap2016_2019.mat','ap')
brewer=[157,183,185]
read_config_

for i=1:3
    ap{i}=[];
    % rs{i}=[];
    for ano=2016:2019
        
        [ap_]=read_ap_obs(brewer(i),strrep(fullfile(Cal.path_root),'2019',num2str(ano)));
        ap{i}=[ap{i};ap_];
        %rs{i}=[rs{i};rs_];
    end
    save('ap2016_2019.mat','-append','ap')


load('ap2016_2019')

ap_ev{i}=meanperiods(ap{i}(:,[1,7:end]), events{i}); 
if i==1
label_ap_1={'T1','T2','T3','HV','V15','V5','V_15','V24','RT','RH','F1','F2','V_5','V8','SLI','SLV'}
label_ap=['Date','year','jday','hour','min','sec',cellfun(@(x) strcat(x,'_off'),label_ap_1,'UniformOutput',false),...
          cellfun(@(x) strcat(x,'_on'),label_ap_1,'UniformOutput',false)];
ap_table{i}=array2table(ap{i},'VariableNames',label_ap');
data=[ap_ev{i}.m(:,2:end),ap_ev{i}.std(:,2:end)]; 
%idx=[1:32,;32+(1:32)];
 l=size(data,2)/2;
 idx=[1:l,;l+(1:l)]
idx=idx(:);
data=data(:,idx);
lbl=[label_ap(7:end),cellfun(@(x) strcat(x,'_std'),label_ap(7:end),'UniformOutput',false)];
lbl=lbl(idx)
ap_ev_table{i}=array2table([ap_ev{i}.m(:,1),ap_ev{i}.N(:,2),data],'VariableNames',['Date','N',lbl]);

else
   label_ap_2={'T1','T2','T3','HV','V12','V5','V_12','V24','RT','T4','T5','T6','V5ss','V_5ss','SLI','SLV',...
           'HGI','HGV','Ex1','Ex2','RH','Hgm3','Ex4','Ex5'}
   label_ap=['Date','year','jday','hour','min','sec',cellfun(@(x) strcat(x,'_off'),label_ap_2,'UniformOutput',false),...
          cellfun(@(x) strcat(x,'_on'),label_ap_2,'UniformOutput',false)];
   ap_table{i}=array2table(ap{i},'VariableNames',label_ap');
   data=[ap_ev{i}.m(:,2:end),ap_ev{i}.std(:,2:end)]; 
   l=size(data,2)/2;
   idx=[1:l,;l+(1:l)];idx=idx(:);
   data=data(:,idx);
   lbl=[label_ap(7:end),cellfun(@(x) strcat(x,'_std'),label_ap(7:end),'UniformOutput',false)];
   lbl=lbl(idx)
   ap_ev_table{i}=array2table([ap_ev{i}.m(:,1),ap_ev{i}.N(:,2),data],'VariableNames',['Date','N',lbl]);

end


 writetable(ap_ev_table{i},'IzoTriad_2016_2019.xls','Sheet',strrep('ap_ev_157','157',num2str(brewer(i))),'WriteRowNames',true)
 writetable(ap_table{i},'IzoTriad_2016_2019.xls','Sheet',strrep('ap_157','157',num2str(brewer(i))),'WriteRowNames',true)
 
 %% Plots for the report
 f(i)=figure(i)
 set(f(i),'Tag','hv_2016_2019');
 ap_avg=ap_ev_table{i};
 apt=ap_table{i};
 
 plot(apt.Date,apt.HV_on,'x',apt.Date,apt.HV_off,':o')
 hold on
 errorbar(ap_avg.Date,ap_avg.HV_on,ap_avg.HV_on_std,'b.')
 datetick
 vline_v(events{i}.dates,'-',strrep(events{i}.labels,'_',' '))
 title(num2str(brewer(i)))
% grid on;
%%
temp=strmatch('T',apt.Properties.VariableNames);
figure
plot(apt.Date,apt{:,temp})
legend(apt.Properties.VariableNames(temp))
clickableLegend(apt.Properties.VariableNames(temp))
datetick('x',12,'keepticks','keeplimits')
grid
vline_v(events{i}.dates,'-',strrep(events{i}.labels,'_',' '))
title(num2str(brewer(i)))
% grid on; 
%% Voltajes
figure
temp=strmatch('V',apt.Properties.VariableNames);
plot(apt.Date,apt{:,temp}-round(mean(apt{:,temp})))
clickableLegend(apt.Properties.VariableNames(temp))
datetick('x',12,'keepticks','keeplimits')
grid
vline_v(events{i}.dates,'-',strrep(events{i}.labels,'_',' '))
ylabel('Voltages difference to the median');
title(num2str(brewer(i)))
% grid on;
% hold on;
% box on;
% errorbar(dt_avg(:,1),dt_avg(:,30),dt_avg(:,31),'r.')
% %h=ploty(smoothdata(dt_avg(:,[1,30,32]),'rloess',90));
% dt_avg_m=meanmonth(dt_avg(:,[1,30,32]));
% h=ploty(dt_avg_m.media(:,[1,5,6]));
% set(h,'Linewidth',4);
% set(h(2),'Color', [0 1 0]);
% set(h(1),'Color', [0 0 0]);
% axis 'tight'
% datetick('x',12,'keepticks','keeplimits')
% ylim(round(median(dt_avg(:,32)))+[-20,20])
% title(num2str(brewer(i)))
% grid on;
% legend('dt low','dt high')    
%     
% end
% 
% for i=1:3
%     
% f(3+i)=figure(3+i)
% set(f(3+i),'Tag','dt_2016_2019');
% % Blacklist filtering
% dt_avg=blacklist_summary(fullfile(Cal.path_root,'..','configs',strcat('blacklist_',num2str(brewer(i)),'.txt')),dt{i});
% % Time interval seleccion
% dateini=datenum('2019-01-01');
% dateend=datenum('2020-01-01');
% dt_avg=dt_avg(dt_avg(:,1) > dateini & dt_avg(:,1) < dateend,:);
% errorbar(dt_avg(:,1),dt_avg(:,32),dt_avg(:,33),'b.')
% hold on;
% box on;
% errorbar(dt_avg(:,1),dt_avg(:,30),dt_avg(:,31),'r.')
% %h=ploty(smoothdata(dt_avg(:,[1,30,32]),'rloess',90));
% dt_avg_m=meanmonth(dt_avg(:,[1,30,32]));
% h=ploty(dt_avg_m.media(:,[1,5,6]));
% set(h,'Linewidth',4);
% set(h(2),'Color', [0 1 0]);
% set(h(1),'Color', [0 0 0]);
% axis 'tight'
% datetick('x',12,'keepticks','keeplimits')
% ylim(round(median(dt_avg(:,32)))+[-20,20])
% title(num2str(brewer(i)))
% grid on;
% legend('dt low','dt high')
end

Width=20; Height=10;
%printfiles_report(f,Cal.dir_figs,'Width',Width,'Height',Height);