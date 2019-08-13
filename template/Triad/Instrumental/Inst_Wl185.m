% options_pub=struct('outputDir',fullfile('..','html'),'showCode',true); 
% close all; publish(fullfile(pwd,'Inst_Wl185.m'),options_pub);

%% Brewer Evaluation
clear all; 
file_setup='calizo2018_setup'; 
run(fullfile('..','..',file_setup));     % configuracion por defecto

Cal.n_inst=find(Cal.brw==185); 
Cal.Date.CALC_DAYS=datenum(2017,6,1):now;

OP_config =Cal.brw_config_files{Cal.n_inst,1};
ALT_config=Cal.brw_config_files{Cal.n_inst,2};

Cal.file_latex = fullfile('..','Triad','latex');
Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures');
mkdir(Cal.dir_figs);

%% Configs: Operative
[~, b c]=fileparts(OP_config);
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
matrix2latex_config(events_cfg_op.data(2:end,:),fullfile(Cal.file_latex,['Op_config_',Cal.brw_str{Cal.n_inst},'.tex']),...
                    'rowlabels',events_cfg_op.legend(2:end),'columnlabels',cellstr(datestr(events_cfg_op.data(1,:),1))',...
                    'size','footnotesize');
                
%%
[c,ia]=unique(events_cfg_op.data(2:end,:)','rows','stable');               
tableop=array2table(c,'VariableNames',str2name(events_cfg_op.legend(2:end)));                
tableop.Date=datetime(datestr(events_cfg_op.data(1,ia)));                
tableop=tableop(:,[end,1:end-1])

%% Configs: Check
[a b c]=fileparts(ALT_config);
try
   if ~strcmpi(strcat(b,c),sprintf('config%d_a.cfg',Cal.brw(Cal.n_inst)))
      fprintf(strcat(1,'\rCUIDADO!! Puede que las configuraciones cargadas no sean las esperadas\n',...
                       '(Expected: %s, Loading: %s)\n'),...
              sprintf('config%d_a.cfg',Cal.brw(Cal.n_inst)),strcat(b,c)); 
   end    
   events_cfg_chk=getcfgs(Cal.Date.CALC_DAYS,ALT_config);    
   fprintf('\nBrewer %s: Alternative Config.\n',Cal.brw_name{Cal.n_inst});
   displaytable(events_cfg_chk.data(2:end,:),cellstr(datestr(events_cfg_chk.data(1,:),1))',12,'.5g',events_cfg_chk.legend(2:end));
catch exception
   fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{Cal.n_inst});
end
matrix2latex_config(events_cfg_chk.data(2:end,:),fullfile(Cal.file_latex,['Chk_config_',Cal.brw_str{Cal.n_inst},'.tex']),...
                    'rowlabels',events_cfg_chk.legend(2:end),'columnlabels',cellstr(datestr(events_cfg_chk.data(1,:),1))',...
                    'size','footnotesize');
                
 %%
[c,ia]=unique(events_cfg_chk.data(2:end,:)','rows','stable')               
tablec=array2table(c,'VariableNames',str2name(events_cfg_op.legend(2:end)))                
tablec.Date=datetime(datestr(events_cfg_op.data(1,ia)));                
tablec=tablec(:,[end,1:end-1])                
                
%%  EVENTOS               
event_info=getevents(Cal,'grp','events');

matrix2latex_ctable(cellstr(datestr(event_info.dates)),fullfile(Cal.file_latex,['table_events_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                   'rowlabels',event_info.labels','columnlabels',{'DATE'},...
                                   'alignment','c','resize',0.9);


event_info.labels
cellstr(datestr(event_info.dates))'

displaytable(event_info.labels',{'events'},30,'%f',cellstr(datestr(event_info.dates))');


%displaytable(cellstr(datestr(event_final.dates)),{'DATE'},14,'.4f',event_final.labels)

%Events importantes
%eventos importantes
ev=table();
ev.label=event_info.labels';
ev.dates=datetime(datestr(event_info.dates));
events_table=ev([1,3,4,5,end],:)

%definidos 19/04/2018
event_info.labels=event_info.labels([1:end]);
event_info.dates=event_info.dates([1:end]);


%event_info.labels=event_info.labels([1,3,5,end-1:end]);
%event_info.dates=event_info.dates([1,3,5,end-1:end]);

% events_table=ev([3,7,8,10],:)
% 
% event_info.labels=event_info.labels([3,7,8,10]);
% event_info.dates=event_info.dates([3,7,8,10]);


% ev185=table();
% ev185.label=event_info.labels';
% ev185.dates=datetime(datestr(event_info.dates));
% r=[4,6,11,13];
% ev185(r,:)=[]
% 
% event_final.labels=event_info.labels;
% event_final.dates=event_info.dates;
% 
% event_final.labels(r)=[];
% event_final.dates(r)=[];
% 
% 
% 
% matrix2latex_ctable(cellstr(datestr(event_final.dates)),fullfile(Cal.file_latex,['table_events_fin_',Cal.brw_str{Cal.n_inst},'.tex']),...
%                                    'rowlabels',event_final.labels','columnlabels',{'DATE'},...
%                                    'alignment','c','resize',0.9);


%% Historical review AVG info
close all; 
%[tabla_avg sl_data]=report_avg(Cal,'grp','events','outlier_flag',{'','','','','','',''});
[tabla_avg,sl_data,dt_data,rs_data,ap_data ]=report_avg(Cal,'grp_custom',event_info,'outlier_flag',{'','','','','','',''});

figure(findobj('Tag','SLAVG_R6')); hold all;
s(1)=stairs([events_cfg_op.data(1,:) sl_data(end,1)],[events_cfg_op.data(17,:) events_cfg_op.data(17,end)],'-','color','m','LineWidth',2);
s(2)=stairs([events_cfg_chk.data(1,:) sl_data(end,1)],[events_cfg_chk.data(17,:) events_cfg_chk.data(17,end)],'-','color','k','LineWidth',2);
%legend(s,{'R6 ref (Op)','R6 ref (Chk)'},'Location','NorthEast');
legend(s,{'R6 ref (Op)'},'Location','NorthEast');
%%
figure(findobj('Tag','SLAVG_R6'));
vline_v(event_info.dates,'k',event_info.labels)

figure(findobj('Tag','SLAVG_F5'));
vline_v(event_info.dates,'k',event_info.labels)
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
set(gca,'Ylim',[25,35])


%%
figure
plot(ap_data(:,1),ap_data(:,4),'-')
vline_v(event_info.dates,'k',event_info.labels);
hold all
plot(ap_data(:,1),ap_data(:,4),'-o')
set(gca,'Ylim',[1000,1150])

datetick

%%
figure
mdt=meanmonth(dt_data);
hold on
shadedErrorBar(mdt.media(:,1),mdt.media(:,8),mdt.sigma(:,8),'r',1)
shadedErrorBar(mdt.media(:,1),mdt.media(:,7),mdt.sigma(:,7),'b',1)
datetick('x',12,'Keeplimits','Keepticks');
hold all;
s(1)=stairs([events_cfg_op.data(1,:) dt_data(end,1)],1E9*[events_cfg_op.data(9,:) events_cfg_op.data(9,end)],'-','color','m','LineWidth',2);
hold on
s(2)=stairs([events_cfg_chk.data(1,:) dt_data(end,1)],1E9*[events_cfg_chk.data(9,:) events_cfg_chk.data(9,end)],'-','color','k','LineWidth',2);
legend(s,{'DT ref (Op)','DT ref (Chk)'},'Location','NorthEast');

set(gcf,'Tag','DT_AVG');
title(['DT  Brewer#', Cal.brw_str{Cal.n_inst}]);
legend(s,{'DT ref (Op)','DT ref (Chk)'},'Location','NorthEast');
box on;grid on;
ylabel('DT ns');
%vline_v(event_info.dates,'k',event_info.labels)



% Command-Print 
fprintf('\r\nAVG data\r\n'); 
if size(tabla_avg.data,1)<15      
   display_table(tabla_avg.data(:,2:end)',tabla_avg.events,15,'.4f',tabla_avg.data_lbl);
else
   display_table(tabla_avg.data(:,2:end),tabla_avg.data_lbl,7,'.4f',tabla_avg.events);
end
matrix2latex_ctable(tabla_avg.data(:,2:end)',fullfile(Cal.file_latex,['table_avg_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                   'rowlabels',tabla_avg.data_lbl,'columnlabels',tabla_avg.events,...
                                   'alignment','c','resize',1.25); %landscape
                               
                               
avg185_table=array2table(tabla_avg.data(:,2:end),'VariableNames',varname(tabla_avg.data_lbl),'RowNames',varname(tabla_avg.events));
avg185_table.Date=datetime(datestr(tabla_avg.data(:,1)))


%% Temperature dependence 
close all; 
t.dates=[datenum(2018,1,1),datenum(2018,2,20),datenum(2018,5,18),datenum(2018,7,18),datenum(2018,8,15)];
t.labels=[{'Start','K&Z maintenance','Power supply','Before Arosa','After Arosa'}];
% [tabla_tc,sl_raw_185]=report_temperature(Cal,OP_config,ALT_config,'grp_custom',t,'reprocess',0);

%[tabla_tc,sl_raw_185]=report_temperature(Cal,OP_config,ALT_config,'grp','month+events');
[tabla_tc,sl_raw_185]=report_temperature(Cal,OP_config,ALT_config,'grp_custom',event_info);
fprintf('\r\nTemp. Dependence\r\n'); 
 if size(tabla_tc.data,1)<9      
    display_table(tabla_tc.data(:,2:end)',tabla_tc.events,12,'.2f',tabla_tc.data_lbl);
 else
    display_table(tabla_tc.data(:,2:end),tabla_tc.data_lbl,8,'.2f',tabla_tc.events);
 end
 matrix2latex_ctable(tabla_tc.data(:,2:end)',fullfile(Cal.file_latex,['table_tc_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_tc.data_lbl,'columnlabels',tabla_tc.events,...
                                    'alignment','c','resize',0.9)
ix=sort(findobj('-regexp','Tag','TEMP_COMP\w+'));
arrayfun( @(x) set(x,'tag',[Cal.brw_str{Cal.n_inst},'_',get(x,'tag')]),ix)
Width=15; Height=20;
printfiles_report(ix',Cal.dir_figs,'Width',Width,'Height',Height);
%close all
%nanmean(tabla_tc.sl{1})
%nanmean(tabla_tc.sl_r{1})

% % %Eventos personalizados: 
% events=struct('dates',[datenum(2016,1,1),datenum(2016,4,16),datenum(2016,5,1)],'labels',{{'PTB','New Lamp','Lamp Stable'}}); 
% %[tabla_tc, sl_rw]=report_temperature(Cal,OP_config,'grp');
% %[tabla_tc]=report_temperature(Cal,OP_config,ALT_config,'grp_custom',event_final);
% [tabla_tc,sl_rw_185]=report_temperature(Cal,OP_config,ALT_config,'grp_custom',events);
% 
% %  Command-Print
% %  %tabla_tc=tabla_alt;
% %  fprintf('\r\nTemp. Dependence\r\n'); 
%  if size(tabla_tc.data,1)<9      
%     displaytable(tabla_tc.data(:,2:end)',tabla_tc.events,12,'.2f',tabla_tc.data_lbl);
%  else
%     displaytable(tabla_tc.data(:,2:end),tabla_tc.data_lbl,8,'.2f',tabla_tc.events);
%  end
% %  matrix2latex_ctable(tabla_tc.data(:,2:end)',fullfile(Cal.file_latex,['table_tc_',Cal.brw_str{Cal.n_inst},'.tex']),...
% %                                     'rowlabels',tabla_tc.data_lbl,'columnlabels',tabla_tc.events,...
% %                                     'alignment','c','resize',0.9);


%% DARK

figure
gscatter(sl_raw_185(:,1),sl_raw_185(:,12),fix(sl_raw_185(:,4)/2.5)*2.5)
datetick
set(gca,'Ylim',[0100,300])
vline_v(event_info.dates,'k',event_info.labels)
title('Dark Count (SL) 185')
grid
%%
figure
hx=gscatter(sl_raw_185(:,4),sl_raw_185(:,12),datestr(fix(sl_raw_185(:,1)/30)*30,'yymm'))
set(gca,'Ylim',[0,300])
title('Dark Count (SL) #185')
grid
xlabel('Temperature')
set(findobj(hx,'DisplayName','1701'),'Marker','x')
set(findobj(hx,'DisplayName','1601'),'Marker','x')


%%

% sl_data_r=cell2mat(tabla_tc.sl_r');
% sl_data=cell2mat(tabla_tc.sl');
% r6=meanperiods(sl_data,events);
% r6_r=meanperiods(sl_data_r,events);
% 
% figure
% mdt=meanmonth(sl_data_r);
% msl=meanmonth(sl_data);
% hold on
% shadedErrorBar(mdt.media(:,1),mdt.media(:,12),mdt.sigma(:,12),'r',1)
% shadedErrorBar(msl.media(:,1),msl.media(:,12),msl.sigma(:,12),'b',1)
% datetick('x',12,'Keeplimits','Keepticks');
% hold all;
% s(1)=stairs([events_cfg_op.data(1,:) sl_data(end,1)],[events_cfg_op.data(17,:) events_cfg_op.data(17,end)],'-','color','m','LineWidth',2);
% hold on
% s(2)=stairs([events_cfg_chk.data(1,:) sl_data(end,1)],[events_cfg_chk.data(17,:) events_cfg_chk.data(17,end)],'-','color','k','LineWidth',2);
% legend(s,{'SL ref (Op)','SL ref (Chk)'},'Location','NorthEast');
% 
% %errorbar(events.dates,r6.m(:,end),r6.std(:,end))
% 
% 
% set(gcf,'Tag','SL_R6');
% title(['SL R6 Brewer#', Cal.brw_str{Cal.n_inst}]);
% legend(s,{'SL ref (Op)','SL ref (Chk)'},'Location','NorthEast');
% box on;grid on;
% ylabel('SL (log counts/sec)');


                                
%%
 %% Temperature dependence alternative configuration
% close all; 
% %[tabla_tc, sl_rw ,Fraw0,Fopt,Forig]=report_temperature(Cal,OP_config,'grp','events');
% %[tabla_tc, sl_rw ,Fraw0,Falt,Forig]=report_temperature(Cal,ALT_config,'grp','events');
% tabla_tc=report_temperature(Cal,ALT_config,'grp');
% % Command-Print
% %tabla_tc=tabla_alt;
%  fprintf('\r\nTemp. Dependence\r\n'); 
%  if size(tabla_tc.data,1)<9      
%     displaytable(tabla_tc.data(:,2:end)',tabla_tc.events,12,'.2f',tabla_tc.data_lbl);
%  else
%     displaytable(tabla_tc.data(:,2:end),tabla_tc.data_lbl,8,'.2f',tabla_tc.events);
%  end
%  matrix2latex_ctable(tabla_tc.data(:,2:end)',fullfile(Cal.file_latex,['table_tc_',Cal.brw_str{Cal.n_inst},'.tex']),...
%                                     'rowlabels',tabla_tc.data_lbl,'columnlabels',tabla_tc.events,...
%                                     'alignment','c','resize',0.9);
                                
%% Filter attenuation
% old filters
close all; 
events_filter=struct('dates',[datenum(2016,1,1),datenum(2017,3,1),datenum(2018,2,4)],'labels',{{'ND original','ND change #1',' ND change #2'}}); 
tabla_fi_old=report_filter(Cal,'grp_custom',events_filter,'date_range',[datenum(2016,1,1),now]);
display_table(tabla_fi_old.data(:,2:end),tabla_fi_old.data_lbl,10,'.2f',tabla_fi_old.events)
%
% Command-Print 
% fprintf('\r\nFI Analysis\r\n'); 
tabla_fi=tabla_fi_old;
 if size(tabla_fi.data,1)<11      
    display_table(tabla_fi.data(:,2:end)',tabla_fi.events,12,'.2f',tabla_fi.data_lbl);
 else
    display_table(tabla_fi.data(:,2:end),tabla_fi.data_lbl,10,'.2f',tabla_fi.events)
 end
% nice graf only for one period
tabla_fi_new=report_filter(Cal,'grp_custom',events_filter,'date_range',[datenum(2018,2,10) now])
%display_table(tabla_fi_new.data(:,2:end),tabla_fi_new.data_lbl,10,'.2f',tabla_fi_new.events)
% Command-Print 
% fprintf('\r\nFI Analysis\r\n');  
% matrix2latex_ctable(tabla_fi.data(:,2:end)',fullfile(Cal.file_latex,['table_fi_',Cal.brw_str{Cal.n_inst},'.tex']),...
%                                     'rowlabels',tabla_fi.data_lbl,'columnlabels',tabla_fi.events,...
%                                     'alignment','c','resize',0.9);
%new filters 
%tabla_fi=report_filter(Cal,'grp','month','date_range',[datenum(2017,3,9):now])% Command-Print 
% fprintf('\r\nFI Analysis\r\n');  
tabla_fi=tabla_fi_new; 
% if size(tabla_fi.data,1)<11      
%     display_table(tabla_fi.data(:,2:end)',tabla_fi.events,12,'.2f',tabla_fi.data_lbl);
%  else
%     display_table(tabla_fi.data(:,2:end),tabla_fi.data_lbl,10,'.2f',tabla_fi.events);
%  end
 matrix2latex_ctable(tabla_fi.data(:,2:end)',fullfile(Cal.file_latex,['table_fi_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                     'rowlabels',tabla_fi.data_lbl,'columnlabels',tabla_fi.events,...
                                     'alignment','c','resize',0.9);

%% SC test

close all; tabla_sc=report_sc(Cal,OP_config,'grp','month+events'); 

% Command-Print 
% fprintf('\r\nSC Analysis\r\n'); 
% if size(tabla_sc.data,1)<10      
%    displaytable(tabla_sc.data(:,2:end)',tabla_sc.events,15,'.0f',tabla_sc.data_lbl);
% else
%    displaytable(tabla_sc.data(:,2:end),tabla_sc.data_lbl,15,'.0f',tabla_sc.events);
% end
% matrix2latex_ctable(tabla_sc.data(:,2:end)',fullfile(Cal.file_latex,['table_sc_',Cal.brw_str{Cal.n_inst},'.tex']),...
%                                     'rowlabels',tabla_sc.data_lbl,'columnlabels',tabla_sc.events,...
%                                     'alignment','c','resize',0.9);

figure; set(gcf,'Tag',strcat('CSN_hist',Cal.brw_str{Cal.n_inst}));
patch([min(tabla_sc.data(:,1))-10 max(tabla_sc.data(:,1))+10 max(tabla_sc.data(:,1))+10 min(tabla_sc.data(:,1))-10],...
      [repmat(1019,1,2) repmat(1021,1,2)], ...
      [.953,.953,.953],'LineStyle',':'); hold on;
errorbar(tabla_sc.data(:,1),tabla_sc.data(:,2),matadd(tabla_sc.data(:,2),-tabla_sc.data(:,3)),matadd(tabla_sc.data(:,4),-tabla_sc.data(:,2)),'*'); 
set(gca,'Layer','Top');
grid; datetick('x',12,'keepLimits','keepTicks');
title(sprintf('Brewer %s (red v line means Hg replacement)',Cal.brw_name{Cal.n_inst})); 
ylabel('Cal Step Number'); box on; axis('tight');
hold on; s=stairs(tabla_sc.data(:,1),tabla_sc.data(:,end-1),'-','color','m','LineWidth',2);
vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
set(vl,'LineStyle','None'); set(findobj(gcf,'Type','Text'),'FontSize',7);
%%
figure;
cs_ref=unique(tabla_sc.data(~isnan(tabla_sc.data(:,5)),5));

set(gcf,'Tag',strcat('CSN_hist',Cal.brw_str{Cal.n_inst}));
patch([min(tabla_sc.data(:,1))-10 max(tabla_sc.data(:,1))+10 max(tabla_sc.data(:,1))+10 min(tabla_sc.data(:,1))-10],...
      [repmat(cs_ref-2,1,2) repmat(cs_ref+2,1,2)], ...
      [.953,.953,.953],'LineStyle',':'); 
datetick('x',12,'keepLimits','keepTicks'); 
hold on;
errorbar(tabla_sc.data(:,1),tabla_sc.data(:,2),matadd(tabla_sc.data(:,2),-tabla_sc.data(:,3)),matadd(tabla_sc.data(:,4),-tabla_sc.data(:,2)),'*'); 
set(gca,'Layer','Top');
grid;
datetick('x',12,'keepLimits','keepTicks');
title(sprintf('Brewer %s (red v line means Hg replacement)',Cal.brw_name{Cal.n_inst})); 
set(gca,'Ylim',cs_ref+[-5,5]);
ylabel('Cal Step Number');
box on; 
%axis('tight');
hold on; 
s=stairs(tabla_sc.data(:,1),tabla_sc.data(:,end-1),'-','color','m','LineWidth',2);
vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
set(vl,'LineStyle','None');
set(findobj(gcf,'Type','Text'),'FontSize',7);


%% CZ Report
close all; [tabla_HS tabla_HL]=report_wavelength(Cal,'grp','events');

% Command-Print
 fprintf('\r\nHS Scans\r\n'); 
 if size(tabla_HS.data,1)<11      
   displaytable(tabla_HS.data(:,2:end)',tabla_HS.events,12,'.4f',tabla_HS.data_lbl);
 else
   displaytable(tabla_HS.data(:,2:end),tabla_HS.data_lbl,11,'.4f',tabla_HS.events);
 end
 matrix2latex_ctable(tabla_HS.data(:,2:end)',fullfile(Cal.file_latex,['table_hs_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_HS.data_lbl,'columnlabels',tabla_HS.events,...
                                    'alignment','c','resize',0.9);
 
 fprintf('\r\nHL Scans\r\n'); 
 if size(tabla_HL.data,1)<11     
   displaytable(tabla_HL.data(:,2:end)',tabla_HL.events,12,'.4f',tabla_HL.data_lbl);
 else
   displaytable(tabla_HL.data(:,2:end),tabla_HL.data_lbl,11,'.4f',tabla_HL.events);
 end
 matrix2latex_ctable(tabla_HL.data(:,2:end)',fullfile(Cal.file_latex,['table_hl_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_HL.data_lbl,'columnlabels',tabla_HL.events,...
                                    'alignment','c','resize',0.9);

%% DSP's 

close all; 
tabla_dsp=report_dispersion(Cal,'grp','month+events','fpath',fullfile(Cal.path_root,'..','DSP'),...
                                           'date_range',datenum(2017,1,1):now,'process',0);  

   events_cfg_op =getcfgs(tabla_dsp.data(:,1)',OP_config);    
   events_cfg_chk=getcfgs(tabla_dsp.data(:,1)',ALT_config);    

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
m185=nanmean(tabla_dsp.data(:,[15,17]))
s185=nanstd(tabla_dsp.data(:,[15,17]))

plot(tabla_dsp.data(:,1),(mq(:,1)-m185(1))/0.001)
hold on
plot(tabla_dsp.data(:,1),(mc(:,1)-m185(2))/0.001)
legend({'cuadratic','cubic'})
plot(tabla_dsp.data(:,1),(tabla_dsp.data(:,15)-m185(1))/0.001,'k+')
hold on
plot(tabla_dsp.data(:,1),(tabla_dsp.data(:,17)-m185(2))/0.001,'bx')

datetick('x','mmm/yy','keepticks')
grid on
ylabel('Micromenter step')
title({'Ozone Abs coeff diff vs mean',sprintf('%s quad=%.4f cubic=%.4f',Cal.brw_str{Cal.n_inst},m185)})
box on
set(gca,'YLim',[-2 2])
set(gca,'XLim',datenum(Cal.Date.cal_year,1,[-180,365]))

vline_v(event_info.dates,'k',event_info.labels);
hfill([-0.5,0.5],'grey')


% %% DSP2
% close all; 
% tabla_dsp=report_dispersion(Cal,'grp','month+events','fpath',fullfile(Cal.path_root,'..','DSP'),...
%                                             'date_range',datenum(2016,1,1):now);  
% % 
% jnan=isnan(tabla_dsp.data(:,2));
% if any(jnan)
% tabla_dsp.data(jnan,:)=[];
% tabla_dsp.events(jnan)=[];
% end
% dsp_table=array2table(tabla_dsp.data(:,2:end),'VariableNames',str2name(tabla_dsp.data_lbl),'Rownames',tabla_dsp.events);
% dsp_table.Date=datetime(datestr(tabla_dsp.data(:,1)));
% 
% 
% 
% 
%    events_cfg_op =getcfgs(tabla_dsp.data(:,1)',OP_config);    
%    events_cfg_chk=getcfgs(tabla_dsp.data(:,1)',ALT_config);    
% % 
% % % Command-Print 
%  fprintf('\r\nDSP Analysis\r\n');  
%  if size(tabla_dsp.data,1)<13      
%     displaytable(tabla_dsp.data(:,[2:8 15:end])',tabla_dsp.events,12,'.4f',tabla_dsp.data_lbl([1:7 14:end]));
%  else
%     displaytable(tabla_dsp.data(:,[2:8 15:end]),tabla_dsp.data_lbl([1:7 14:end]),10,'.4f',tabla_dsp.events);
%  end
%  matrix2latex_ctable(tabla_dsp.data(:,15:end),fullfile(Cal.file_latex,['table_dsp_',Cal.brw_str{Cal.n_inst},'.tex']),...
%                                     'rowlabels',tabla_dsp.events,'columnlabels',tabla_dsp.data_lbl(14:end),...
%                                     'alignment','c','resize',0.9);
% % 
% figure; set(gcf,'Tag',strcat('DSP_A1',Cal.brw_str{Cal.n_inst})); 
% tabla_dsp.data(isnan(tabla_dsp.data(:,2)),:)=[];
% l1=ploty(tabla_dsp.data(:,[1 15 17]),'s');
% datetick('x',12,'KeepLimits','keepTicks'); 
% lim=get(gca,'XLim'); set(gca,'YLim',[0.3365 0.3455],'XLim',[lim(1)-10 lim(2)+10]); 
% set(l1(1),'Marker','s','LineStyle','-','Color','k','MarkerFaceColor','k'); 
% set(l1(2),'Marker','s','LineStyle','-','Color','g','MarkerFaceColor','g');
% 
% hold on; s(2)=stairs([events_cfg_chk.data(1,:) tabla_dsp.data(end,1)],[events_cfg_chk.data(7,:) events_cfg_chk.data(7,end)],'-','color','m','LineWidth',2);
%          s(1)=stairs([events_cfg_op.data(1,:) tabla_dsp.data(end,1)],[events_cfg_op.data(7,:) events_cfg_op.data(7,end)],'-','color','b','LineWidth',2);
% %vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
% vline_v(event_final.dates,'',event_final.labels)
% %set(vl,'LineStyle','None'); set(findobj(gca,'Type','text'),'FontSize',7);    
% title(sprintf('A1: %s',Cal.brw_name{Cal.n_inst}));
% 
% ylabel('Ozone Absorption Coefficient (A1)'); grid
% set(findobj(gcf,'Type','text'),'BackgroundColor','w');
% 
% m185=nanmean(tabla_dsp.data(:,[15,17]));
% hline(m185(1),'k',sprintf('quad=%.4f ',m185(1)))
% hline(m185(2),'g',sprintf(' cubic=%.4f',m185(2)))
% legendflex([l1;s(1);s(2)],{'A1 Quad','A1 Cubic','A1 ref (Op)','A1 ref (Chk)'}, ...
%        'anchor', {'sw','sw'}, 'buffer',[15 15], 'nrow',2,'fontsize',12,'box','off');
% 
% %
% figure
% mq=mean_smooth_abs(tabla_dsp.data(:,1),tabla_dsp.data(:,[15]),90,0);
% mc=mean_smooth_abs(tabla_dsp.data(:,1),tabla_dsp.data(:,[17]),90,0);
% hold on
% m185=nanmean(tabla_dsp.data(:,[15,17]));
% 
% 
% plot(tabla_dsp.data(:,1),(mq(:,1)-m185(1))/0.001)
% hold on
% plot(tabla_dsp.data(:,1),(mc(:,1)-m185(2))/0.001)
% legend({'cuadratic','cubic'})
% plot(tabla_dsp.data(:,1),(tabla_dsp.data(:,15)-m185(1))/0.001,'k+')
% hold on
% plot(tabla_dsp.data(:,1),(tabla_dsp.data(:,17)-m185(2))/0.001,'bx')
% 
% datetick('x','mmm/yy','keepticks')
% grid on
% ylabel('Micromenter step')
% title({'Brewer  #185 Ozone Abs coeff diff vs mean',sprintf('quad=%.4f cubic=%.4f',m185)})
% box on
% set(gca,'YLim',[-2 2])
% vline_v(event_info.dates,'k',event_info.labels);
% hfill([-0.5,0.5],'grey')
% % printfiles_report(gcf,Cal.dir_figs,'Width',18,'Height',12,'LineMode','Scaled');