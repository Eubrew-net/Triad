%% GENERA GRAFICOS Y TABLAS PARA EL REPORT
%  tablas de temperatura
%  
% read obs from Triad
sl=cell(3,1);sl_r=cell(3,1);sl_s=cell(3,1);sl_raw=cell(3,1);
dsum=cell(3,1);
dsum_r=cell(3,1);
brewer=[157,183,185]

read_config_
%load('events_16_19.mat')
close all
%%
% eventos temperatura 157
ev157=[    {'09-Mar-2016'}    {'17-Mar-2016'}    {'22-May-2017'}    {'31-Mar-2018'}    {'17-Oct-2018'}    {'10-Dec-2018'}    {'15-Apr-2019'}  {'24-May-2019'}];
t{1}.dates=datenum(ev157)';
t{1}.labels=[{'Kipp_Zonen'}    {'Cal1000W'}    {'SLJump'}    {'DT'}    {'roof157'}    {'SLReplacement'}    {'SLJump_1'}    {'FixPowerSuply'}];
%eventos temperatura 183
t{2}.dates=[datenum(2015,6,1),datenum(2017,12,10),datenum(2018,05,12)];
t{2}.labels=[{'start','hv','SL chg'}];
%eventos de temperatura 185
t{3}.dates=[datenum(2015,6,10),datenum(2015,10,15),datenum(2015,11,30),datenum(2016,1,16),datenum(2016,2,1),datenum(2016,3,10),datenum(2016,4,10),datenum(2018,2,14),datenum(2018,3,26),datenum(2019,7,1)];
t{3}.labels={'IOS','sl_jump','SC','PTB','IZO','KZ','SLnew','HV','HVS','Huelva'};
%% 
%  Analizamos los datos
for i=1:3
 Cal.n_inst=i;   
 [tabla_tc,sl_raw_157]=report_temperature(Cal,OP_config,ALT_config,'grp_custom',t{i},'reprocess',0);
 %% salidas latex
  matrix2latex_ctable(tabla_tc.data(:,2:end)',fullfile(Cal.file_latex,['table_tc_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_tc.data_lbl,'columnlabels',tabla_tc.events,...
                                    'alignment','c','resize',0.9);
                             
   ix=sort(findobj('-regexp','Tag','TEMP_COMP\w+'));
   arrayfun( @(x) set(x,'tag',[Cal.brw_str{Cal.n_inst},'_',get(x,'tag')]),ix)
   %Width=20; Height=15;
   printfiles_report(ix',Cal.dir_figs,'Width',20,'Height',15);
 
end

%% Figuras de temperatura
%
for i=1:3
    sl{i}=[];
    sl_s{i}=[];
    sl_r{i}=[];
    sl{i}.tabla=[];
    sl{i}.events={};
    sl_raw{i}=[];
    dsum_r{i}=[];
    dsum{i}=[];
    
    s1_=strrep(fullfile(Cal.path_root,'..','Triad','Instrumental','IZO#157_sl_rw.mat'),'157',num2str(brewer(i)))
    if exist(s1_)
        s=load(s1_);
        sl_raw{i}=[sl_raw{i};s.sl_rw];
        sl{i}.tabla=[sl{i}.tabla;s.tabla_tc.data]
        %sl{i}.data=[sl{i}.data;s.tabla_tc.data];
        sl_s{i}=[sl_s{i},s.tabla_tc.sl];
        sl_r{i}=[sl_s{i},s.tabla_tc.sl_r];
        sl{i}.labels=[sl{i}.events,s.tabla_tc.events];
        sl{i}.dates=cellfun(@(x) min(x(:,1)),s.tabla_tc.sl_r);
    end    
    hold off;
    figure;
    for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[1,end]),'.'), end; grid; title(num2str(brewer(i)));
    axis tight; datetick('keeplimits')
    vline_v(sl{i}.dates,'-',sl{i}.labels)
    grid on; box on;
    title(Cal.brw_str(i))
    
for ano=2016:2019
    j=ano-2015;
       
     s1_=(strrep( strrep(fullfile(Cal.path_root,'Triad','Langley','summary_old_Brw157_2019.txt'),'2019',num2str(ano)),'157',num2str(brewer(i))))
     s2_=(strrep( strrep(fullfile(Cal.path_root,'Triad','Langley','summary_Brw157_2019.txt'),'2019',num2str(ano)),'157',num2str(brewer(i))))
 
     ds=load(s1_);
     ds_r=load(s2_);
     dsum_r{i}=[dsum_r{i};ds_r(:,1:15)]; % some summaries include config
     dsum{i}=[dsum{i};ds(:,1:15)];
     
  
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
    %writetable(sl_ev{i},'sl_2016_2019.xls','Sheet',num2str(brewer(i)),'WriteRowNames',true)
   
   
   %% temperature analysis periods
   
    

   f1=figure
   h1=plot(dsum{i}(:,1),dsum{i}(:,end),'r+',dsum{i}(:,1),dsum{i}(:,end-1),'k-');
   hold on
   h2=plot(dsum_r{i}(:,1),dsum_r{i}(:,end),'mo',dsum_r{i}(:,1),dsum_r{i}(:,end-1),'b-');
   vline_v(events{i}.dates,'-',strrep(events{i}.labels,'_',' '))
   legend([h1;h2],{'R6 op','op_ref','R6 chk','chk_ref'})
   set(h1(2),'LineWidth',2);
   set(h2(2),'LineWidth',2);
   
   grid;
   title(num2str(brewer(i)))
   datetick('x',12,'keepticks','keeplimits')
   xlim([datenum(2015,12,1),now])
   set(gca,'LooseInset',get(gca,'TightInset'))
    
   
   
   % operative vs temp
   f2=figure;
   for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[2,end]),'.'); end; grid; title(num2str(brewer(i)));
   %  alternative vs temp
   f3=figure;
   for j=1:length(sl_r{i}), hold all, ploty(sl_r{i}{j}(:,[2,end]),'.'); end; grid; title(num2str(brewer(i)));

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

%% load configuration
%read_config_
%load alt_cfg_16_19
%load op_cfg_16_19
%events{i}=op_cfg{1}.Properties.RowNames;
%read_config_
%load('events_16_19.mat')

%events{i}=op_cfg{1}.Properties.RowNames;
