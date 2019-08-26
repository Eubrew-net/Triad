% options_pub=struct('outputDir',fullfile('..','html'),'showCode',true); 
% close all; publish(fullfile(pwd,'Inst_Wl157_2016_2019.m'),options_pub);

%% Brewer Evaluation
clear all; 
file_setup='calizo_setup'; 
run(fullfile('..',file_setup));     % configuracion por defecto

Cal.n_inst=find(Cal.brw==157); 
Cal.Date.CALC_DAYS=datenum(2015,6,1):now;

OP_config =Cal.brw_config_files{Cal.n_inst,1};
ALT_config=Cal.brw_config_files{Cal.n_inst,2};

Cal.file_latex = fullfile('..','latex');


Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures');
mkdir(Cal.dir_figs);

Cal.dir_tables   = fullfile(Cal.file_latex,filesep(),'tables');
mkdir(Cal.dir_tables);


%% Configs: Operative
[a b c]=fileparts(OP_config);
try
   if ~strcmpi(strcat(b,c),sprintf('config%d.cfg',Cal.brw(Cal.n_inst)))
      fprintf(strcat(1,'\rCUIDADO!! Puede que las configuraciones cargadas no sean las esperadas\n',...
                       '(Expected: %s, Loading: %s)\n'),...
              sprintf('config%d.cfg',Cal.brw(Cal.n_inst)),strcat(b,c)); 
   end
   events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,OP_config);    
   fprintf('\nBrewer %s: Operative Config.\n',Cal.brw_name{Cal.n_inst});
   displaytable(events_cfg_op.data(2:end,:),cellstr(datestr(events_cfg_op.data(1,:),1))',12,'.5g',events_cfg_op.legend(2:end));
catch exception
   fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{Cal.n_inst});
end
matrix2latex_config(events_cfg_op.data(2:end,:),fullfile(Cal.dir_tables,['Op_config_',Cal.brw_str{Cal.n_inst},'.tex']),...
                    'rowlabels',events_cfg_op.legend(2:end),'columnlabels',cellstr(datestr(events_cfg_op.data(1,:),1))',...
                    'size','footnotesize');
                
[c,ia]=unique(events_cfg_op.data(2:end,:)','rows','stable');               
tableop_157=array2table(c,'VariableNames',str2name(events_cfg_op.legend(2:end)));                
tableop_157.Date=datetime(datestr(events_cfg_op.data(1,ia)));                
tableop_157=tableop_157(:,[end,1:end-1])
%% Configs: Check
[a b c]=fileparts(ALT_config);
try
   if ~strcmpi(strcat(b,c),sprintf('config%d_a.cfg',Cal.brw(Cal.n_inst)))
      fprintf(strcat(1,'\rCUIDADO!! Puede que las configuraciones cargadas no sean las esperadas\n',...
                       '(Expected: %s, Loading: %s)\n'),...
              sprintf('config%d_a.cfg',Cal.brw(Cal.n_inst)),strcat(b,c)); 
   end    
   events_cfg_chk=getcfgs(Cal.Date.CALC_DAYS,ALT_config);    
   fprintf('\nBrewer %s: Second Config.\n',Cal.brw_name{Cal.n_inst});
   %displaytable(events_cfg_chk.data(2:end,:),cellstr(datestr(events_cfg_chk.data(1,:),1))',12,'.5g',events_cfg_chk.legend(2:end));
   t157=display_table(events_cfg_chk.data(1:end,:),cellstr(datestr(events_cfg_chk.data(1,:),1))',12,'.5g',events_cfg_chk.legend(1:end));
   
catch exception
   fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{Cal.n_inst});
end
matrix2latex_config(events_cfg_chk.data(2:end,:),fullfile(Cal.file_latex,['Chk_config_',Cal.brw_str{Cal.n_inst},'.tex']),...
                    'rowlabels',events_cfg_chk.legend(2:end),'columnlabels',cellstr(datestr(events_cfg_chk.data(1,:),1))',...
                    'size','footnotesize');
                
[c,ia]=unique(events_cfg_chk.data(2:end,:)','rows','stable')               
tablec_157=array2table(c,'VariableNames',str2name(events_cfg_op.legend(2:end)))                
tablec_157.Date=datetime(datestr(events_cfg_op.data(1,ia)));                
tablec_157=tablec_157(:,[end,1:end-1])                
                
%%
%%  EVENTOS               
event_info=getevents(Cal,'gr','events');


matrix2latex_ctable(cellstr(datestr(event_info.dates)),fullfile(Cal.file_latex,['table_events_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                   'rowlabels',event_info.labels','columnlabels',{'DATE'},...
                                   'alignment','c','resize',0.9);


event_info.labels
cellstr(datestr(event_info.dates))'

t_events=display_table(event_info.labels',{'events'},30,'%f',cellstr(datestr(event_info.dates))');
t_events.Date=event_info.dates;
t_events.Fecha=datetime(datestr(event_info.dates))

%%eventos importantes-> SELECCIONAMOS LOS IMPORTANTES<-
ev=table();
ev.label=event_info.labels';
ev.dates=datetime(datestr(event_info.dates))
events_table=ev;


event_info.labels=event_info.labels;
event_info.dates=event_info.dates;


% event_info.labels=event_info.labels(4:end);
% event_info.dates=event_info.dates(4:end);
 displaytable(event_info.labels',{'events'},30,'%f',cellstr(datestr(event_info.dates))');                


                
               

%% Historical review AVG info
close all; 
[tabla_avg,sl_data,dt_data,rs_data,ap_data ]=report_avg(Cal,'grp','events','outlier_flag',{'sl','','','','','',''});      

f=figure(findobj('Tag','SLAVG_F5'));
vline_v(event_info.dates,'k',event_info.labels)
hold all;

f=figure(findobj('Tag','SLAVG_R6'));
vline_v(event_info.dates,'k',event_info.labels)
 set(gca,'YLim',median(events_cfg_op.data(17,:))+[-30,30])
hold all;
s(1)=stairs([events_cfg_op.data(1,:) sl_data(end,1)],[events_cfg_op.data(17,:) events_cfg_op.data(17,end)],'-','color','m','LineWidth',2);
hold on
s(2)=stairs([events_cfg_chk.data(1,:) sl_data(end,1)],[events_cfg_chk.data(17,:) events_cfg_chk.data(17,end)],'-','color','k','LineWidth',2);
legend(s,{'R6 ref (Op)','R6 ref (Chk)'},'Location','NorthEast');

f=figure(maxf(findobj('tag','SLAVG_R6')));
Width=36; Height=12;
R6tag=[Cal.brw_str{Cal.n_inst},'_',get(f,'tag')];
set(f,'tag',R6tag)
ylim([320 400])

f=figure(maxf(findobj('tag','DTAVG')));
hold on
try
 [a,b]=findchangepts(rmmissing(dt_data(:,5)),'MinThreshold', 1700,'MinDistance',7,'Statistic','mean');
 
 vline_v(datenum(dt_data(a,1)),'r',cellstr(datestr(dt_data(a,1))))
catch
    disp('signal');
end


printfiles_report(f,Cal.dir_figs,'Width',Width,'Height',Height);
                           
%%
table_avg_157=array2table(tabla_avg.data,'VariableNames',str2name(['fecha',tabla_avg.data_lbl]));               
table_avg_157.Date=datetime(datestr(tabla_avg.data(:,1)));   
table_avg_157.Properties.RowNames=varname(tabla_avg.events);   
table_avg_157=table_avg_157(:,[end,1:end-1])      
try
  rows2vars(table_avg_157)
catch
  disp('transpose');
end

input.data=table_avg_157(:,2:end); %(solo los nï¿½meros,quitamos las fechas()
input.transposeTable = 1;
input.tableCaption = 'Brewer #156 average table';
latex_avg = latexTable(input);
cellwrite('table_avg_157.tex',latex_avg);
writetable(table_avg_157,'table_avg_157.xls')


%%
 figure(findobj('Tag','SLAVG_F5'));
 vline_v(event_info.dates,'k',event_info.labels)
 
%%

figure(findobj('-regexp','Tag',R6tag));
vline_v(event_info.dates,'k',event_info.labels)

try
 [a,b]=findchangepts(rmmissing(sl_data(:,12)),'MinThreshold', 50,'MinDistance',7,'Statistic','mean') 
 vline_v(datenum(sl_data(a,1)),'r',cellstr(datestr(sl_data(a,1))))
catch
    disp('need signal');
end 

%%
  figure
  plot(ap_data(:,1),ap_data(:,6),'x')
  vline_v(event_info.dates,'r',event_info.labels)
  datetick('x','yy/mmm')
  ylabel('Amps');
  title('SL lamp Intesity');
%% DT
figure
mdt=meanmonth(dt_data,12);
hold on
shadedErrorBar(mdt.media(:,1),mdt.media(:,8),mdt.sigma(:,8),'r',1)
shadedErrorBar(mdt.media(:,1),mdt.media(:,7),mdt.sigma(:,7),'b',1)
datetick('x',12,'Keeplimits','Keepticks');
gx=get(gca,'XLim');
hold all;

%DT table

  dt_op=[[events_cfg_op.data(1,:) dt_data(end,1)];1E9*[events_cfg_op.data(9,:) events_cfg_op.data(9,end)]]';
  dt_alt=[[events_cfg_chk.data(1,:) dt_data(end,1)];1E9*[events_cfg_chk.data(9,:) events_cfg_chk.data(9,end)]]';
  

s(1)=stairs([events_cfg_op.data(1,:) dt_data(end,1)],1E9*[events_cfg_op.data(9,:) events_cfg_op.data(9,end)],'-','color','m','LineWidth',2);
hold on
s(2)=stairs([events_cfg_chk.data(1,:) dt_data(end,1)],1E9*[events_cfg_chk.data(9,:) events_cfg_chk.data(9,end)],'-','color','k','LineWidth',2);
legend(s,{'DT ref (Op)','DT ref (Chk)'},'Location','NorthEast');

set(gcf,'Tag','DT_AVG');
title(['DT  Brewer#', Cal.brw_str{Cal.n_inst}]);
legend(s,{'DT ref (Op)','DT ref (Chk)'},'Location','NorthEast');
box on
ylabel('DT ns');
set(gca,'Xlim',gx)


%%
%%
figure(maxf(findobj('Tag','DTAVG')));
hold all;
s=[];
vline_v(event_info.dates,'k',event_info.labels)
s(1)=stairs(tabla_avg.data(:,1),tabla_avg.data(:,7),'linewidth',4);
s(2)=stairs(tabla_avg.data(:,1),tabla_avg.data(:,9),'linewidth',4);
s(3)=stairs([events_cfg_op.data(1,:) sl_data(end,1)],1E9*[events_cfg_op.data(9,:) events_cfg_op.data(9,end)],'-','color','m','LineWidth',2);
s(4)=stairs([events_cfg_chk.data(1,:) sl_data(end,1)],1E9*([events_cfg_chk.data(9,:) events_cfg_chk.data(9,end)]),'-','color','k','LineWidth',2);
legend(s,{'DT high','Dt low','DT ref (Op)','DT ref (Chk)'},'Location','NorthEast');
%set(gca,'Ylim',[25,35])
f=figure(maxf(findobj('tag','DTAVG')));
Width=15; Height=20;
set(f,'tag',[Cal.brw_str{Cal.n_inst},'_',get(f,'tag')])
printfiles_report(f,Cal.dir_figs,'Width',Width,'Height',Height);


%% Command-Print 
fprintf('\r\nAVG data\r\n'); 
% if size(tabla_avg.data,1)<15      
%    displaytable(tabla_avg.data(:,2:end)',tabla_avg.events,15,'.4f',tabla_avg.data_lbl);
% else
%    displaytable(tabla_avg.data(:,2:end),tabla_avg.data_lbl,7,'.4f',tabla_avg.events);
% end
avg_t=display_table(tabla_avg.data(:,2:end)',tabla_avg.events,15,'.4f',tabla_avg.data_lbl);
writetable(avg_t,'avg_ev_2019_157.txt');
matrix2latex_ctable(tabla_avg.data(:,2:end)',fullfile(Cal.dir_tables,['table_avg_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                   'rowlabels',tabla_avg.data_lbl,'columnlabels',str2var(tabla_avg.events),...
                                   'alignment','c','resize',0.9);
                               
                               
 avg_t=array2table(tabla_avg.data(:,1:end)','RowNames',['Date',str2name(tabla_avg.data_lbl)],'VariableNames',varname(cellstr(datestr(tabla_avg.data(:,1)))))  ;                           
 %avg_t.Fecha=datetime(datestr(avg_t.Date));
 %avg_t=[avg_t,tabla_avg.events];

%input.data =avg_t(:,2:11);
% lx=latexTable(input);
% lx=strrep(lx,'_',' ');
% lx=strrep(lx,'x ',' ');
% fid=fopen('sl_avg_157.tex','w');
% [nrows,ncols] = size(lx);
% for row = 1:nrows
%     fprintf(fid,'%s\n',lx{row,:});
% end
% fclose(fid);

%avg_t.Date=datetime(datestr(tabla_avg.data(:,1)))

%% Temperature dependence

close all; 
 
t.dates=[datenum(2015,6,1),datenum(2016,02,01),datenum(2016,04,10),datenum(2018,2,15),datenum(2018,3,25)];
t.labels=[{'ios','x1','kzm','j1','j2'}];
%vline_v(t.dates,'r-',t.labels)
%t.dates=[datenum(2019,1,1)];
%t.labels=[{'2019'}];
%  Hay que mirar si hay cambios en la intensidad de la lampara ->
%  
%  aparece una mejora que no se confirma cuando se pone en practica ¿? 
t_temp=t_events([4,5,11,14,17,18,20,23]+2,:)
t.dates=t_temp.Date;
t.labels=t_temp.events;

[tabla_tc,sl_raw_157]=report_temperature(Cal,OP_config,ALT_config,'grp_custom',t,'reprocess',0);
%[tabla_tc,sl_raw_157]=report_temperature(Cal,OP_config,ALT_config,'grp','events','reprocess',0)
%% Global temperature
sl=cat(1,tabla_tc.sl{:});
sl_r=cat(1,tabla_tc.sl_r{:});
sl_m=grpstats(sl(:,[1,end]),{year(sl(:,1)),month(sl(:,1))});
sl_t=grpstats(sl(:,[2,end]),{year(sl(:,1)),month(sl(:,1))});
slr_m=grpstats(sl_r(:,[1,end]),{year(sl(:,1)),month(sl(:,1))});
slr_t=grpstats(sl_r(:,[2,end]),{year(sl(:,1)),month(sl(:,1))});


f=figure;
set(f,'Tag','SL_R6');
plot(sl(:,1),sl(:,end),'.');hold all
plot(sl_m(:,1), sl_m(:,2),'lineWidth',4)


datetick('keeplimits');
sx=get(gca,'XLim');
vline_v(event_info.dates,'k',event_info.labels)
% set(gca,'YLim',median(events_cfg_op.data(17,:))+[-30,30])
s(1)=stairs([events_cfg_op.data(1,:) sl_data(end,1)],[events_cfg_op.data(17,:) events_cfg_op.data(17,end)],'-','color','m','LineWidth',2);
hold on
s(2)=stairs([events_cfg_chk.data(1,:) sl_data(end,1)],[events_cfg_chk.data(17,:) events_cfg_chk.data(17,end)],'-','color','k','LineWidth',2);
set(gca,'XLim',sx);
legend(s,{'obs','mean','R6 ref (Op)','R6 ref (Chk)'});
Width=24; Height=12;
set(f,'tag',[Cal.brw_str{Cal.n_inst},'_',get(f,'tag')])
set(gca,'XLim',sx);
printfiles_report(f,Cal.dir_figs,'Width',Width,'Height',Height,'no_export');


%%
%tabla_tc_157=display_table(tabla_tc.data(:,2:end)',tabla_tc.events,12,'.2f',tabla_tc.data_lbl)

%[tabla_tc,sl_raw_157]=report_temperature(Cal,OP_config,ALT_config,'grp','events');
% Command-Print
%  fprintf('\r\nTemp. Dependence\r\n'); 
%  if size(tabla_tc.data,1)<9      
%     displaytable(tabla_tc.data(:,2:end)',tabla_tc.events,12,'.2f',tabla_tc.data_lbl);
%  else
%     displaytable(tabla_tc.data(:,2:end),tabla_tc.data_lbl,8,'.2f',tabla_tc.events);
%  end 
 matrix2latex_ctable(tabla_tc.data(:,2:end)',fullfile(Cal.file_latex,['table_tc_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_tc.data_lbl,'columnlabels',tabla_tc.events,...
                                    'alignment','c','resize',0.9);
                                
ix=sort(findobj('-regexp','Tag','TEMP_COMP\w+'));
arrayfun( @(x) set(x,'tag',[Cal.brw_str{Cal.n_inst},'_',get(x,'tag')]),ix)

printfiles_report(ix',Cal.dir_figs,'Width',Width,'Height',Height);
close all

nanmean(tabla_tc.sl{1})
nanmean(tabla_tc.sl_r{1})

%% SL reference 371
% 
% 
% config_temp.n_inst=1;
% config_temp.brw_name='#157';
% config_temp.final_days=datenum(2016,1,1):now;
% 
% NTC={}; tabla_regress={}; ajuste={}; Args={};
% 
% s=load('IZO#157_sl_rw.mat')
% 
% %% Temperature dependence. 
%     
% %     [NTC{1},ajuste{1},Args{1},Fr]=temp_coeff_raw(config_temp,s.sl_rw,'outlier_flag',1,'intensity_flag',1);
% %     
% %     disp(sprintf(' ORIG MS9: %5.0f +/-%2.0f  %3.1f +/- %3.2f  ',ajuste{1}.orig(7,[1 3 2 4])));
% %     disp(sprintf('  NEW MS9: %5.0f +/-%2.0f  %3.1f +/- %3.2f  ',ajuste{1}.new(7,[1 3 2 4])));
%     dtime=[t.dates(1),t.dates(2)];
%     TCorig=unique(tableop_157{2:end,2:6},'rows')
%     [NTC1,ajustex,Argsx,Fraw,Forig]=temp_coeff_raw(config_temp,s.sl_rw,'outlier_flag',1,'plots',1,...
%                                 'N_TC',TCorig(1:5)','date_range',dtime);
%     NTC{1}=tabla_tc.data(1,end-4:end)                       
%     [NTC2,ajustex,Argsx,Fraw,Fnew]=temp_coeff_raw(config_temp,s.sl_rw,'outlier_flag',1,'plots',0,...
%                                 'N_TC',NTC{1},'date_range',dtime);
% Forigx=Forig; Fn=Fnew;
% 
% figure;  set(gcf,'Tag','TEMP_COMP_DATE');
% plot(Forigx(:,1),Forigx(:,2),'b.','MarkerSize',6); 
% ylabel('Temperature','Color','b'); ax(1)=gca; set(ax(1),'YAxisLocation','right','XTicklabels',{' '}); 
% [mn,sn]=grpstats(Forigx(:,[1,end,2]),{year(Forigx(:,1)),fix(Forigx(:,1))},{'mean','sem'});
% [mt,st]=grpstats(Fn(:,[1,end]),{year(Fn(:,1)),fix(Fn(:,1))},{'mean','sem'}); 
% ax(2) = axes('YAxisLocation','left','Color','none'); 
% hold all; errorbar(mn(:,1),mn(:,2),sn(:,2),'Color','k','Marker','s');
% errorbar(mt(:,1),mt(:,2),st(:,2),'Color','g','Marker','s');
% errorbar(mt(:,1),mt(:,2),st(:,2),'Color','g','Marker','s');
% title(['R6 Temperature dependence Brewer#', Cal.brw_str{Cal.n_inst}]); ylabel('Standard Lamp R6 ratio');
% datetick('x',6,'KeepTicks','KeepLimits'); grid on; 
% lg=legend(ax(2),'Old temperature coeff','New temperature coeff','Location','best'); 
% set(lg,'HandleVisibility','Off');  set(findobj(lg,'Type','text'),'FontSize',7,'HandleVisibility','Off');    
% linkprop(ax,{'Position','XTick'}); 
% 
% figure; set(gcf,'Tag','TEMP_COMP_TEMP')
% [mn,sn]=grpstats(Forigx(:,[2,end]),Forigx(:,2),{'mean','sem'});
% [mt,st]=grpstats(Fn(:,[2,end]),Fn(:,2),{'mean','sem'});
% hold on; errorbar(mn(:,1),mn(:,2),sn(:,2),'Color','k','Marker','s');
% errorbar(mt(:,1),mt(:,2),st(:,2),'Color','g','Marker','s'); grid;
% lg=legend('Old temperature coeff','New temperature coeff','Location','best');
% set(lg,'HandleVisibility','Off'); set(findobj(lg,'Type','text'),'FontSize',7,'HandleVisibility','Off');    
% title(['R6 Temperature dependence Brewer#', Cal.brw_str{Cal.n_inst}]);
% ylabel('Standard Lamp R6 ratio'); xlabel('Temperature'); set(gca,'Box','On');


 
%% dark evaluation
figure
ploty(sl_raw_157(:,[4,12]),'+')
title('Dark Count (SL) 157')
xlabel('Temperature')
grid

%% DARK

figure
gscatter(sl_raw_157(:,1),sl_raw_157(:,12),fix(sl_raw_157(:,4)/2.5)*2.5)
datetick
set(gca,'Ylim',[0,60])
vline_v(event_info.dates,'k',event_info.labels)
title('Dark Count (SL) 157')
grid
%%
figure
hx=gscatter(sl_raw_157(:,4),sl_raw_157(:,12),datestr(fix(sl_raw_157(:,1)/30)*30,'yymm'))
set(gca,'Ylim',[000,60])
title('Dark Count (SL) #157')
grid
xlabel('Temperature')
% set(findobj(hx,'DisplayName','1701'),'Marker','x')
% set(findobj(hx,'DisplayName','1601'),'Marker','x')


                                
%% Filter attenuation

% para evauar por periodos, por ejemplos cuando se cambia un filtro
close all; 

t.dates=[datenum(2010,6,1),datenum(2018,6,1),datenum(2019,4,12)];
t.labels=[{'historic','last year','now'}];

tabla_fi=report_filter(Cal,'grp_custom',t,'date_range',t.dates);


% Command-Print 
% fprintf('\r\nFI Analysis\r\n');  
%  if size(tabla_fi.data,1)<11      
%     displaytable(tabla_fi.data(:,2:end)',tabla_fi.events,12,'.2f',tabla_fi.data_lbl);
%  else
%     displaytable(tabla_fi.data(:,2:end),tabla_fi.data_lbl,10,'.2f',tabla_fi.events);
%  end
 displaytable(tabla_fi.data(:,2:end),tabla_fi.data_lbl,10,'.2f',tabla_fi.events);
 matrix2latex_ctable(tabla_fi.data(:,2:end)',fullfile(Cal.file_latex,['table_fi_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                     'rowlabels',str2latex(tabla_fi.data_lbl),'columnlabels',str2latex(tabla_fi.events),...
                                     'alignment','c','resize',0.9);

%% SC's 
close all; 

Cal.Date.CALC_DAYS=datenum(2015,6,1):now;
tabla_sc=report_sc(Cal,OP_config,'grp','events');
tabla_sc_157=display_table(tabla_sc.data(:,2:end),tabla_sc.data_lbl,15,'.0f',tabla_sc.events)
writetable(tabla_sc_157,'SC_summary_2019_157.txt')
writetable(tabla_sc_157,'SC_summary_2019_157.xls')

% Command-Print 
% fprintf('\r\nSC Analysis\r\n'); 
% if size(tabla_sc.data,1)<10      
%    displaytable(tabla_sc.data(:,2:end)',tabla_sc.events,15,'.0f',tabla_sc.data_lbl);
% else
%    displaytable(tabla_sc.data(:,2:end),tabla_sc.data_lbl,15,'.0f',tabla_sc.events);
% end
 matrix2latex_ctable(tabla_sc.data(:,2:end)',fullfile(Cal.file_latex,['table_sc_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                     'rowlabels',str2latex(tabla_sc.data_lbl),'columnlabels',str2latex(tabla_sc.events),...
                                     'alignment','c','resize',0.9);

                                
figure; set(gcf,'Tag',strcat('SC_CSN_hist',Cal.brw_str{Cal.n_inst}));
patch([min(tabla_sc.data(:,1))-10 max(tabla_sc.data(:,1))+10 max(tabla_sc.data(:,1))+10 min(tabla_sc.data(:,1))-10],...
      [repmat(1025,1,2) repmat(1027,1,2)], ...
      [.953,.953,.953],'LineStyle',':'); hold on;
errorbar(tabla_sc.data(:,1),tabla_sc.data(:,2),matadd(tabla_sc.data(:,2),-tabla_sc.data(:,3)),matadd(tabla_sc.data(:,4),-tabla_sc.data(:,2)),'*'); 
set(gca,'Layer','Top');
grid; datetick('x',12,'keepLimits','keepTicks');
title(sprintf('Brewer %s ' ,Cal.brw_name{Cal.n_inst})); 
ylabel('Cal Step Number'); box on; axis('tight');
hold on; s=stairs(tabla_sc.data(:,1),tabla_sc.data(:,end-1),'-','color','m','LineWidth',2);
vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
set(vl,'LineStyle','None'); set(findobj(gcf,'Type','Text'),'FontSize',7);

ix=sort(findobj('-regexp','Tag','SC\w+'));
arrayfun( @(x) set(x,'tag',[Cal.brw_str{Cal.n_inst},'_',get(x,'tag')]),ix)

printfiles_report(ix',Cal.dir_figs,'Width',Width,'Height',Height);
close all
writetable(tabla_sc_157,'SC_summary_157.txt')
writetable(tabla_sc_157,'SC_summary_157.xls')

%nanmean(tabla_tc.sl{1})
%nanmean(tabla_tc.sl_r{1})



%% CZ Report
close all; [tabla_HS tabla_HL]=report_wavelength(Cal,'grp','events');
tabla_hl_157=display_table(tabla_HS.data(:,2:end)',tabla_HS.events,12,'.4f',tabla_HS.data_lbl)
tabla_hs_157=display_table(tabla_HL.data(:,2:end)',tabla_HL.events,12,'.4f',tabla_HL.data_lbl)

% Command-Print
%  fprintf('\r\nHS Scans\r\n'); 
%  if size(tabla_HS.data,1)<11      
%    displaytable(tabla_HS.data(:,2:end)',tabla_HS.events,12,'.4f',tabla_HS.data_lbl);
%  else
%    displaytable(tabla_HS.data(:,2:end),tabla_HS.data_lbl,11,'.4f',tabla_HS.events);
%  end
  matrix2latex_ctable(tabla_HS.data(:,2:end)',fullfile(Cal.file_latex,['table_hs_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                     'rowlabels',tabla_HS.data_lbl,'columnlabels',tabla_HS.events,...
                                     'alignment','c','resize',0.9);
 
%  fprintf('\r\nHL Scans\r\n'); 
%  if size(tabla_HL.data,1)<11     
%    displaytable(tabla_HL.data(:,2:end)',tabla_HL.events,12,'.4f',tabla_HL.data_lbl);
%  else
%    displaytable(tabla_HL.data(:,2:end),tabla_HL.data_lbl,11,'.4f',tabla_HL.events);
%  end
 matrix2latex_ctable(tabla_HL.data(:,2:end)',fullfile(Cal.file_latex,['table_hl_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_HL.data_lbl,'columnlabels',tabla_HL.events,...
                                    'alignment','c','resize',0.9);

%% DSP's 
close all; 
Cal.Date.CALC_DAYS=datenum(2015,6,1):now;
[tabla_dsp,dsp_quad,dsp_cubic]=report_dispersion_new(Cal,'grp','events','fpath',fullfile(Cal.path_root,'..','DSP'),'process',0,...
                                           'date_range',Cal.Date.CALC_DAYS);  

Cal.Date.CALC_DAYS=datenum(2015,6,1):now;
[tabla_dsp_week,dsp_quad,dsp_cubic]=report_dispersion_new(Cal,'grp','week','fpath',fullfile(Cal.path_root,'..','DSP'),'process',0,...
                                           'date_range',Cal.Date.CALC_DAYS);  
                                       
fprintf('\r\nDSP  week Analysis\r\n');  
jn=~isnan(tabla_dsp_week.data(:,2));
dsp_week=array2table(tabla_dsp_week.data(jn,1:end),'VariableNames',varname(tabla_dsp_week.data_lbl),'Rownames',varname(tabla_dsp_week.events(jn)));
dsp_week.Date=datetime(datestr(tabla_dsp_week.data(jn,1)));

writetable(dsp_week,'dsp_2019_157.txt')


t_dsp_157=display_table(tabla_dsp_week.data(jn,[1,2 15,17]),tabla_dsp_week.data_lbl([1,2 15,17]),10,'.4f',tabla_dsp_week.events(jn));
matrix2latex_ctable(tabla_dsp_week.data(jn,15:end),fullfile(Cal.file_latex,['table_dsp_week',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_dsp_week.events,'columnlabels',tabla_dsp_week.data_lbl(14:end),...
                                    'alignment','c','resize',0.9);
                                      
                                       
                                       

events_cfg_op =getcfgs(tabla_dsp.data(:,1)',OP_config);    
events_cfg_chk=getcfgs(tabla_dsp.data(:,1)',ALT_config);    
jn=~isnan(tabla_dsp.data(:,2));
%tabla_dsp_157=display_table(tabla_dsp.data(jn,[2:7 14:end]),tabla_dsp.data_lbl([2:7 14:end]),10,'.4f',tabla_dsp.events(jn))

dsp_ev=array2table(tabla_dsp.data(jn,1:end),'VariableNames',varname(tabla_dsp.data_lbl),'Rownames',varname(tabla_dsp.events(jn)));
dsp_ev.Date=datetime(datestr(tabla_dsp.data(jn,1)));

writetable(dsp_ev,'dsp_avg_2019_157.txt')


% Command-Print 
fprintf('\r\nDSP Analysis\r\n'); 
jn=~isnan(tabla_dsp.data(:,2));
if size(tabla_dsp.data,1)<13      
    displaytable(tabla_dsp.data(jn,[2:8 15:end])',tabla_dsp.events(jn),12,'.4f',tabla_dsp.data_lbl([2:7 14:end]));
 else
    displaytable(tabla_dsp.data(:,[2:8 15:end]),tabla_dsp.data_lbl([2:7 14:end]),10,'.4f',tabla_dsp.events);
 end
 matrix2latex_ctable(tabla_dsp.data(jn,15:end),fullfile(Cal.file_latex,['table_dsp_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',str2latex(tabla_dsp.events(jn)),'columnlabels',str2latex(tabla_dsp.data_lbl(14:end)),...
                                    'alignment','c','resize',0.9);

                                %%   Historic 
figure; set(gcf,'Tag',strcat('DSP_A1_HIST',Cal.brw_str{Cal.n_inst})); 
tabla_dsp.data(isnan(tabla_dsp.data(:,2)),:)=[];
l1=ploty(tabla_dsp.data(:,[1 15 17]),'s');
set(l1(1),'Marker','s','LineStyle','-','Color','k','MarkerFaceColor','k'); 
set(l1(2),'Marker','s','LineStyle','-','Color','g','MarkerFaceColor','g');
lim=get(gca,'XLim'); set(gca,'YLim',[0.3360 0.3435],'XLim',[lim(1)-10 lim(2)+10]); 
hold on; s(2)=stairs([events_cfg_chk.data(1,:) tabla_dsp.data(end,1)],[events_cfg_chk.data(7,:) events_cfg_chk.data(7,end)],'-','color','m','LineWidth',2);
         s(1)=stairs([events_cfg_op.data(1,:) tabla_dsp.data(end,1)],[events_cfg_op.data(7,:) events_cfg_op.data(7,end)],'-','color','b','LineWidth',2);
vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
set(vl,'LineStyle','None'); set(findobj(gca,'Type','text'),'FontSize',7);    
title(sprintf('A1: %s',Cal.brw_name{Cal.n_inst}));
datetick('x',12,'KeepLimits','keepTicks'); ylabel('Ozone Absorption Coefficient (A1)'); grid
legendflex([l1;s(1);s(2)],{'A1 Quad','A1 Cubic','A1 ref (Op)','A1 ref (Chk)'}, ...
                           'anchor', {'sw','sw'}, 'buffer',[10 -5], 'nrow',2,'fontsize',8,'box','on');
set(findobj(gcf,'Type','text'),'BackgroundColor','w');

printfiles_report(gcf,Cal.dir_figs,'Width',21,'Height',9,'LineMode','Scaled');

%%

                            
                                                                     
figure; 
set(gcf,'Tag',strcat('DSP_A1',Cal.brw_str{Cal.n_inst})); 
tabla_dsp.data(isnan(tabla_dsp.data(:,2)),:)=[];
l1=ploty(tabla_dsp.data(:,[1 15 17]),'s');
set(l1(1),'Marker','s','LineStyle','-','Color','k','MarkerFaceColor','k'); 
set(l1(2),'Marker','s','LineStyle','-','Color','g','MarkerFaceColor','g');
lim=get(gca,'XLim');
set(gca,'YLim',[0.3360 0.3435],'XLim',[lim(1)-10 lim(2)+10]); 
hold on; 
s(2)=stairs([events_cfg_chk.data(1,:) tabla_dsp.data(end,1)],[events_cfg_chk.data(7,:) events_cfg_chk.data(7,end)],'-','color','m','LineWidth',2);
 s(1)=stairs([events_cfg_op.data(1,:) tabla_dsp.data(end,1)],[events_cfg_op.data(7,:) events_cfg_op.data(7,end)],'-','color','b','LineWidth',2);

%vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
%set(vl,'LineStyle','None'); set(findobj(gca,'Type','text'),'FontSize',7);    

title(sprintf('A1: %s',Cal.brw_name{Cal.n_inst}));
datetick('x',12,'KeepLimits','keepTicks'); ylabel('Ozone Absorption Coefficient (A1)'); grid
legendflex([l1;s(1);s(2)],{'A1 Quad','A1 Cubic','A1 ref (Op)','A1 ref (Chk)'}, ...
                           'anchor', {'sw','sw'}, 'buffer',[10 -5], 'nrow',2,'fontsize',12,'box','off');
set(findobj(gcf,'Type','text'),'BackgroundColor','w');

printfiles_report(gcf,Cal.dir_figs,'Width',14,'Height',7,'LineMode','Scaled');



%%
fprintf('\r\nDSP Analysis\r\n');  
jn=~isnan(tabla_dsp.data(:,2));
% if any(jnan)
% tabla_dsp.data(jnan,:)=[];
% tabla_~dsp.events(jnan)=[];
% end
dsp_table=array2table(tabla_dsp.data(jn,1:end),'VariableNames',varname(tabla_dsp.data_lbl),'Rownames',varname(tabla_dsp.events(jn)));
dsp_table.Date=datetime(datestr(tabla_dsp.data(:,1)));

t_dsp_157=display_table(tabla_dsp.data(jn,[1,2 15,17]),tabla_dsp.data_lbl([1,2 15,17]),10,'.4f',tabla_dsp.events(jn))


matrix2latex_ctable(tabla_dsp.data(jn,15:end),fullfile(Cal.file_latex,['table_dsp_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_dsp.events,'columnlabels',tabla_dsp.data_lbl(14:end),...
                                    'alignment','c','resize',0.9);
                                
                                

                                      


%%
figure; 
set(gcf,'Tag','DSP_A1_AVG')
l1=ploty(tabla_dsp.data(:,[1 15 17]),'s'); 
hold all
errorbar(tabla_dsp.data(:,1),tabla_dsp.data(:,15),tabla_dsp.data(:,16));
errorbar(tabla_dsp.data(:,1),tabla_dsp.data(:,17),tabla_dsp.data(:,18));
set(l1(1),'Marker','s','LineStyle','-','Color','k','MarkerFaceColor','k'); 
set(l1(2),'Marker','s','LineStyle','-','Color','g','MarkerFaceColor','g');
lim=get(gca,'XLim');
set(gca,'YLim',[0.3360 0.3435],'XLim',[lim(1)-10 lim(2)+10]); 
vline_v(event_info.dates,'r',event_info.labels)
hold on; 
s(2)=stairs([events_cfg_chk.data(1,:) tabla_dsp.data(end,1)],[events_cfg_chk.data(7,:) events_cfg_chk.data(7,end)],'-','color','m','LineWidth',2);
s(1)=stairs([events_cfg_op.data(1,:) tabla_dsp.data(end,1)],[events_cfg_op.data(7,:) events_cfg_op.data(7,end)],'-','color','b','LineWidth',2);

%vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
%set(vl,'LineStyle','None'); set(findobj(gca,'Type','text'),'FontSize',7);    

title(sprintf('A1: %s',Cal.brw_name{Cal.n_inst}));
datetick('x',12,'KeepLimits','keepTicks'); ylabel('Ozone Absorption Coefficient (A1)'); grid
set(findobj(gcf,'Type','text'),'BackgroundColor','w');
m157=nanmean(tabla_dsp.data(:,[15,17]));
hline(m157(1),'k',sprintf('quad=%.4f ',m157(1)))
hline(m157(2),'g',sprintf(' cubic=%.4f',m157(2)))
legendflex([l1;s(1);s(2)],{'A1 Quad','A1 Cubic','A1 ref (Op)','A1 ref (Chk)'}, ...
       'anchor', {'sw','sw'}, 'buffer',[15 15], 'nrow',2,'fontsize',12,'box','off');

%%
figure
set(gcf,'Tag','DSP_A1_RES')
mq=mean_smooth_abs(tabla_dsp.data(:,1),tabla_dsp.data(:,[15]),90,0);
mc=mean_smooth_abs(tabla_dsp.data(:,1),tabla_dsp.data(:,[17]),90,0);
hold on
m157=nanmean(tabla_dsp.data(:,[15,17]))
s157=nanstd(tabla_dsp.data(:,[15,17]))


plot(tabla_dsp.data(:,1),(mq(:,1)-m157(1))/0.001)
hold on
plot(tabla_dsp.data(:,1),(mc(:,1)-m157(2))/0.001)
legend({'cuadratic','cubic'})
plot(tabla_dsp.data(:,1),(tabla_dsp.data(:,15)-m157(1))/0.001,'k+')
hold on
plot(tabla_dsp.data(:,1),(tabla_dsp.data(:,17)-m157(2))/0.001,'bx')

datetick('x','mmm/yy','keepticks')
grid on
ylabel('Micromenter step')
title({'Brewer  #157 Ozone Abs coeff diff vs mean',sprintf('quad=%.4f cubic=%.4f',m157)})
box on
set(gca,'YLim',[-2 2])
vline_v(event_info.dates,'k',event_info.labels);
hfill([-0.5,0.5],'grey')
