clear all
close all
sl=cell(3,1);sl_r=cell(3,1);sl_s=cell(3,1);sl_raw=cell(3,1);
dsum=cell(3,1);
dsum_r=cell(3,1);
brewer=[157,183,185]

read_config_
%load('events_16_19.mat')
close all


for i=1:3
    sl{i}=[];
    sl_s{i}=[];
    sl_r{i}=[];
    sl{i}.tabla=[];
    sl{i}.events={};
    sl_raw{i}=[];
    dsum_r{i}=[];
    dsum{i}=[];
for ano=2019:2019
    j=ano-2018;
    %s1_=fullfile(Cal.path_root(1:end-5),'Triad','instrumental',strcat('IZO#',num2str(brewer(i)),'_sl_rw.mat'))
    s1_=(strrep( strrep(fullfile(Cal.path_root,'Triad','Instrumental','IZO#157_sl_rw.mat'),'2019',num2str(ano)),'157',num2str(brewer(i))))
   if exist(s1_)
     s=load(s1_);
     sl_raw{i}=[sl_raw{i};s.sl_rw];
     sl{i}.tabla=[sl{i}.tabla;s.tabla_tc.data]
     %sl{i}.data=[sl{i}.data;s.tabla_tc.data];
     sl_s{i}=[sl_s{i},s.tabla_tc.sl];
     sl_r{i}=[sl_s{i},s.tabla_tc.sl_r];
     sl{i}.events=[sl{i}.events,s.tabla_tc.events];
     
     %s1_=fullfile(Cal.path_root(1:end-5),num2str(ano),'Triad','Langley',strcat('summary_old_Brw',num2str(brewer(i)),'_',num2str(ano),'.txt'))
     %s2_=fullfile(Cal.path_root(1:end-5),num2str(ano),'Triad','Langley',strcat('summary_Brw',num2str(brewer(i)),'_',num2str(ano),'.txt'))
     s1_=(strrep( strrep(fullfile(Cal.path_root,'Triad','Langley','summary_old_Brw157_2019.txt'),'2019',num2str(ano)),'157',num2str(brewer(i))))
     s2_=(strrep( strrep(fullfile(Cal.path_root,'Triad','Langley','summary_Brw157_2019.txt'),'2019',num2str(ano)),'157',num2str(brewer(i))))
 
     ds=load(s1_);
     ds_r=load(s2_);
     dsum_r{i}=[dsum_r{i};ds_r(:,1:15)]; % some summaries include config
     dsum{i}=[dsum{i};ds(:,1:15)];
     
    
   end
  
end
try
   t_temperature=array2table(sl{i}.tabla,'VariableNames',[{'date'},str2name((s.tabla_tc.data_lbl))],'RowNames',varname(sl{i}.events));
catch
     t_temperature=array2table(sl{i}.tabla,'VariableNames',[{'date'},str2name((s.tabla_tc.data_lbl))]);
     warning('wrong events')
end
   t_temperature.Date_str=datestr(t_temperature.date,'yyyy/mm/dd');
   t_temperature=t_temperature(:,[end,1:end-1])
   writetable(t_temperature,'IzoTriad_2016_2019.xls','Sheet',strrep('temp_157','157',num2str(brewer(i))),'WriteRowNames',true)
   
    sl_o_brw=meanperiods(dsum{i}, events{i}); 
    sl_a_brw=meanperiods(dsum_r{i}, events{i}); 
    data=[events{i}.dates round(sl_o_brw.m(:,end),1) round(sl_o_brw.std(:,end),2)  op_cfg{i}{17,:}' ...
      round(sl_a_brw.m(:,end),1)     round(sl_a_brw.std(:,end),2)  alt_cfg{i}{17,:}' (sl_a_brw.N(:,end))];
    sl_ev{i}=array2table(data,'VariableNames',{'Date','SL','std','SL_op_ref','SL_r','std_r','SL_chk_ref','N'},'RowNames',varname(str2name(strrep(events{i}.labels,'"',''))));
    sl_ev{i}.Fecha=datetime(datestr(sl_ev{i}.Date));
    sl_ev{i}=timetable2table(table2timetable(sl_ev{i}));
    writetable(sl_ev{i},'IzoTriad_2016_2019.xls','Sheet',strrep('sl_157','157',num2str(brewer(i))),'WriteRowNames',true)
    %writetable(sl_ev{i},'sl_2019.xls','Sheet',num2str(brewer(i)),'WriteRowNames',true)
   
   
   %% temperature analysis periods
   
%    hold off;figure(i*4+1);for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[1,end]),'.'), end; grid; title(num2str(brewer(i)));
%    axis tight; datetick('keeplimits')
%    vline_v(sl{i}.tabla(:,1),'-',sl{i}.events)
%    

   f((i-1)*3+1)=figure(i*3+1)
   set(f((i-1)*3+1),'Tag','sl_2019');
   h1=plot(dsum{i}(:,1),dsum{i}(:,end),'r+',dsum{i}(:,1),dsum{i}(:,end-1),'k-');
   hold on
   h2=plot(dsum_r{i}(:,1),dsum_r{i}(:,end),'mo',dsum_r{i}(:,1),dsum_r{i}(:,end-1),'b-');
   vline_v(events{i}.dates,'-',events{i}.labels)
   legend([h1;h2],{'R6 op','op_ref','R6 chk','chk_ref'})
   set(h1(2),'LineWidth',2);
   set(h2(2),'LineWidth',2);
   
   grid;
   title(num2str(brewer(i)))
   datetick('x',12,'keepticks','keeplimits')
   xlim([datenum(2018,12,1),now])
   set(gca,'LooseInset',get(gca,'TightInset'))
    
   
   
   % operative vs temp
   f((i-1)*3+2)=figure(i*3+2);for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[2,end]),'.'); end; grid; title(num2str(brewer(i)));
   set(f((i-1)*3+2),'Tag','sl_2019');
   %  alternative vs temp
   f((i-1)*3+3)=figure(i*3+3);for j=1:length(sl_r{i}), hold all, ploty(sl_r{i}{j}(:,[2,end]),'.'); end; grid; title(num2str(brewer(i)));
   set(f((i-1)*3+3),'Tag','sl_2019');

%    figure(i*4+3);for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[1,end]),'.'), end; grid; title(num2str(brewer(i)))
%     axis tight; datetick('keeplimits')
%    vline_v(sl{i}.tabla(:,1),'-',sl{i}.events)
%
%    figure(i*4+4)
%    plot(dsum{i}(:,1),dsum{i}(:,end),'r-',dsum{i}(:,1),dsum{i}(:,end-1),'k-')
%    hold on
%    plot(dsum_r{i}(:,1),dsum_r{i}(:,end),'m-',dsum_r{i}(:,1),dsum_r{i}(:,end-1),'b.')
%    vline_v(events{i}.dates,'-',events{i}.labels)
%    legend({'R6 op','op_ref','R6 chk','chk_ref'})
%    grid; title(num2str(brewer(i)))
end

Width=20; Height=10;
printfiles_report(f,Cal.dir_figs,'Width',Width,'Height',Height);

%% load configuration
%read_config_
%load alt_cfg_16_19
%load op_cfg_16_19
%events{i}=op_cfg{1}.Properties.RowNames;
%read_config_
%load('events_16_19.mat')

%events{i}=op_cfg{1}.Properties.RowNames;
