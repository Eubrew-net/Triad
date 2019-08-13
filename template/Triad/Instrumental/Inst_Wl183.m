% options_pub=struct('outputDir',fullfile('..','html'),'showCode',true); 
% close all; publish(fullfile(pwd,'Inst_Wl183.m'),options_pub);

%% Brewer Evaluation
clear all; 
ii=2
file_setup='calizo2018_setup'; 
run(fullfile('..','..',file_setup));     % configuracion por defecto

Cal.n_inst=find(Cal.brw==183);
Cal.Date.CALC_DAYS=datenum(2017,6,1):now;

OP_config =Cal.brw_config_files{Cal.n_inst,1};
ALT_config=Cal.brw_config_files{Cal.n_inst,2};

Cal.file_latex = fullfile('latex');
Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures');
mkdir(Cal.dir_figs);

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
matrix2latex_config(events_cfg_op.data(2:end,:),fullfile(Cal.file_latex,['Op_config_',Cal.brw_str{Cal.n_inst},'.tex']),...
                    'rowlabels',events_cfg_op.legend(2:end),'columnlabels',cellstr(datestr(events_cfg_op.data(1,:),1))',...
                    'size','footnotesize');
%%
[c,ia]=unique(events_cfg_op.data(2:end,:)','rows','stable');               
tableop_183=array2table(c,'VariableNames',str2name(events_cfg_op.legend(2:end)));                
tableop_183.Date=datetime(datestr(events_cfg_op.data(1,ia)));                
tableop_183=tableop_183(:,[end,1:end-1])
   

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
   displaytable(events_cfg_chk.data(2:end,:),cellstr(datestr(events_cfg_chk.data(1,:),1))',12,'.5g',events_cfg_chk.legend(2:end));
catch exception
   fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{Cal.n_inst});
end
matrix2latex_config(events_cfg_chk.data(2:end,:),fullfile(Cal.file_latex,['Chk_config_',Cal.brw_str{Cal.n_inst},'.tex']),...
                    'rowlabels',events_cfg_chk.legend(2:end),'columnlabels',cellstr(datestr(events_cfg_chk.data(1,:),1))',...
                    'size','footnotesize');
%% falta etiquetarlos
[c,ia]=unique(events_cfg_chk.data(2:end,:)','rows','stable')               
tablec_183=array2table(c,'VariableNames',str2name(events_cfg_op.legend(2:end)));                
tablec_183.Date=datetime(datestr(events_cfg_op.data(1,ia)));                
tablec_183=tablec_183(:,[end,1:end-1])

%%  EVENTOS               
event_info=getevents(Cal,'grp','events');

matrix2latex_ctable(cellstr(datestr(event_info.dates)),fullfile(Cal.file_latex,['table_events_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                   'rowlabels',event_info.labels','columnlabels',{'DATE'},...
                                   'alignment','c','resize',0.9);


event_info.labels
cellstr(datestr(event_info.dates))'

display_table(event_info.labels',{'events'},30,'%f',cellstr(datestr(event_info.dates))');


%%eventos importantes
% ia: eventos que cambiaron las configuracion
% [c,ia]=unique(events_cfg_chk.data(2:end,:)','rows','stable')               

ev=table();
ev.labels=event_info.labels';
ev.date=datetime(datestr(event_info.dates));
ev.dates=event_info.dates;
events_table=ev([1,2,3,4,end-1:end],:)
%events_table=ev(ia,:)
%events_table=ev(:,:)

event_info=table2struct(events_table,'ToScalar',true)
%event_info.labels=event_info.labels([1,3,4,end-1:end]);
%event_info.dates=event_info.dates([1,3,4,end-1:end]);
%displaytable(event_info.labels',{'events'},30,'%f',cellstr(datestr(event_info.dates))');                

%% Historical review AVG info
close all; 
%[tabla_avg,sl_data,dt_data,rs_data,ap_data ]=report_avg(Cal,'grp_custom',event_info,'outlier_flag',{'sl','dt','','','','',''},'date_range',[datenum(2017,1,1),now]);
%[tabla_avg,sl_data,dt_data,rs_data,ap_data ]=report_avg(Cal,'grp_custom',events,'outlier_flag',{'sl','dt','','','','',''},'date_range',[datenum(2017,1,1),now]);
[tabla_avg,sl_data,dt_data,rs_data,ap_data ]=report_avg(Cal,'grp','events','outlier_flag',{'','dt','','','','',''});      

figure(findobj('Tag','SLAVG_R6')); hold all;
s(1)=stairs([events_cfg_op.data(1,:) sl_data(end,1)],[events_cfg_op.data(17,:) events_cfg_op.data(17,end)],'-','color','m','LineWidth',2);
s(2)=stairs([events_cfg_chk.data(1,:) sl_data(end,1)],[events_cfg_chk.data(17,:) events_cfg_chk.data(17,end)],'-','color','k','LineWidth',2);
legend(s,{'R6 ref (Op)','R6 ref (Chk)'},'Location','NorthEast');

figure(findobj('Tag','SLAVG_R6'));
vline_v(event_info.dates,'k',event_info.labels)

figure(findobj('Tag','SLAVG_F5'));
vline_v(event_info.dates,'k',event_info.labels)

%%
figure(maxf(findobj('Tag','DTAVG')));
set(gca,'Ylim',[19,26])
hold all;
vline_v(event_info.dates,'k',event_info.labels)
s(1)=stairs(tabla_avg.data(:,1),tabla_avg.data(:,7),'linewidth',4);
s(2)=stairs(tabla_avg.data(:,1),tabla_avg.data(:,9),'linewidth',4);
s(3)=stairs([events_cfg_op.data(1,:) sl_data(end,1)],1E9*[events_cfg_op.data(9,:) events_cfg_op.data(9,end)],'-','color','m','LineWidth',2);
s(4)=stairs([events_cfg_chk.data(1,:) sl_data(end,1)],1E9*([events_cfg_chk.data(9,:) events_cfg_chk.data(9,end)]),'-','color','k','LineWidth',2);
legend(s,{'DT high','Dt low','DT ref (Op)','DT ref (Chk)'},'Location','NorthEast');


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
vline_v(event_info.dates,'k',event_info.labels)




%% Voltage events
% que occurre el dia 5 de Noviembre ¿?
figure
plot(ap_data(:,1),ap_data(:,4))
vline_v(event_info.dates,'k',event_info.labels)
datetick
set(gca,'YLim',[1530,1550])
set(gca,'XLim',datenum(2017,10,[1,300]))

h=vline_v(datenum(2017,11,5),'r-',datestr(datenum(2017,11,5))); set(h,'LineWidth',2)
h1=vline_v(datenum(2018,1,14),'r-',datestr(datenum(2018,1,14))); set(h1,'LineWidth',2)
set(gca,'Tag','VoltEvents');
grid
printfiles_report(gcf,Cal.dir_figs,'LineMode','Scaled');

% %%
% figure
% mdt=meanmonth(dt_data);
% hold on
% shadedErrorBar(mdt.media(:,1),mdt.media(:,8),mdt.sigma(:,8),'r',1)
% shadedErrorBar(mdt.media(:,1),mdt.media(:,7),mdt.sigma(:,7),'b',1)
% datetick('x',12,'Keeplimits','Keepticks');
% hold all;
% s(1)=stairs([events_cfg_op.data(1,:) dt_data(end,1)],1E9*[events_cfg_op.data(9,:) events_cfg_op.data(9,end)],'-','color','m','LineWidth',2);
% hold on
% s(2)=stairs([events_cfg_chk.data(1,:) dt_data(end,1)],1E9*[events_cfg_chk.data(9,:) events_cfg_chk.data(9,end)],'-','color','k','LineWidth',2);
% legend(s,{'DT ref (Op)','DT ref (Chk)'},'Location','NorthEast');
% 
% set(gcf,'Tag','DT_AVG');
% title(['DT  Brewer#', Cal.brw_str{Cal.n_inst}]);
% legend(s,{'DT ref (Op)','DT ref (Chk)'},'Location','NorthEast');
% box on;grid on;
% ylabel('DT ns');
% vline_v(event_info.dates,'k',event_info.labels)



%% Command-Print 
fprintf('\r\nAVG data\r\n'); 
if size(tabla_avg.data,1)<15      
   displaytable(tabla_avg.data(:,2:end)',tabla_avg.events,15,'.4f',tabla_avg.data_lbl);
else
   displaytable(tabla_avg.data(:,2:end),tabla_avg.data_lbl,7,'.4f',tabla_avg.events);
end
matrix2latex_ctable(tabla_avg.data(:,2:end)',fullfile(Cal.file_latex,['table_avg_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                   'rowlabels',tabla_avg.data_lbl,'columnlabels',tabla_avg.events,...
                                   'alignment','c','resize',0.9);
 display_table(tabla_avg.data(:,2:end)',tabla_avg.events,15,'.4f',tabla_avg.data_lbl)
%% Temperature dependence
t.dates=[datenum(2018,1,14),datenum(2018,2,1),datenum(2018,4,18)];
t.labels=[{'hv recover','tower','New SL'}];
[tabla_tc_183,sl_raw_183]=report_temperature(Cal,OP_config,ALT_config,'grp_custom',t,'reprocess',0);
%[tabla_tc,sl_raw_183]=report_temperature(Cal,OP_config,ALT_config,'grp','events');

 fprintf('\r\nTemp. Dependence\r\n'); 
 tabla_tc=tabla_tc_183;
 if size(tabla_tc.data,1)<9      
    displaytable(tabla_tc.data(:,2:end)',tabla_tc.events,12,'.2f',tabla_tc.data_lbl);
 else
    displaytable(tabla_tc.data(:,2:end),tabla_tc.data_lbl,8,'.2f',tabla_tc.events);
 end
 matrix2latex_ctable(tabla_tc.data(:,2:end)',fullfile(Cal.file_latex,['table_tc_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_tc.data_lbl,'columnlabels',tabla_tc.events,...
                                    'alignment','c','resize',0.9);

%% Filter attenuation
close all; 
tabla_fi=report_filter(Cal,'grp_custom',event_info,'date_range',[datenum(2015,2,1) now]);
set(gca,'YLim',[-20,10])
set(gca,'XLim',datenum(2016,1,[1,800]))


%tabla_fi=report_filter(Cal,'grp_custom',event_info,'date_range',[datenum(2016,2,1) now]);
% Command-Print 
%fprintf('\r\nFI Analysis\r\n');  
%if size(tabla_fi.data,1)<11      
%   displaytable(tabla_fi.data(:,2:end),tabla_fi.events,12,'.2f',tabla_fi.data_lbl);
%else
%   displaytable(tabla_fi.data(:,2:end),tabla_fi.data_lbl,10,'.2f',tabla_fi.events);
%end
display_table(tabla_fi.data(:,2:end),tabla_fi.data_lbl,10,'.2f',tabla_fi.events)
matrix2latex_ctable(tabla_fi.data(:,2:end)',fullfile(Cal.file_latex,['table_fi_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_fi.data_lbl,'columnlabels',tabla_fi.events,...
                                    'alignment','c','resize',0.9);

%% SC's 
close all; tabla_sc=report_sc(Cal,OP_config,'grp','month+events'); %close all;

% Command-Print 
fprintf('\r\nSC Analysis\r\n'); 
if size(tabla_sc.data,1)<10      
   displaytable(tabla_sc.data(:,2:end)',tabla_sc.events,15,'.0f',tabla_sc.data_lbl);
else
   displaytable(tabla_sc.data(:,2:end),tabla_sc.data_lbl,15,'.0f',tabla_sc.events);
end
matrix2latex_ctable(tabla_sc.data(:,2:end)',fullfile(Cal.file_latex,['table_sc_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_sc.data_lbl,'columnlabels',tabla_sc.events,...
                                    'alignment','c','resize',0.9);
                                
%%                                

figure; set(gcf,'Tag',strcat('CSN_hist',Cal.brw_str{Cal.n_inst}));
patch([min(tabla_sc.data(:,1))-10 max(tabla_sc.data(:,1))+10 max(tabla_sc.data(:,1))+10 min(tabla_sc.data(:,1))-10],...
      [repmat(1019,1,2) repmat(1021,1,2)], ...
      [.953,.953,.953],'LineStyle',':'); hold on;
errorbar(tabla_sc.data(:,1),tabla_sc.data(:,2),matadd(tabla_sc.data(:,2),-tabla_sc.data(:,3)),matadd(tabla_sc.data(:,4),-tabla_sc.data(:,2)),'*'); 
%set(gca,'Layer','Top');
%set(gca,'XLim',datenum(2017,1,[1,350]))


grid; datetick('x',12,'keepLimits','keepTicks');
title(sprintf('Brewer %s (red v line means Hg replacement)',Cal.brw_name{Cal.n_inst})); 
ylabel('Cal Step Number'); box on; axis('tight');
hold on; s=stairs(tabla_sc.data(:,1),tabla_sc.data(:,end-1),'-','color','m','LineWidth',2);
vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
set(vl,'LineStyle','None'); set(findobj(gcf,'Type','Text'),'FontSize',7);
set(gca,'YLim',[nanmean(tabla_sc.data(:,2))-5 nanmean(tabla_sc.data(:,2))+5])
set(gca,'XLim',[min(tabla_sc.data(:,1))-10 nanmax(tabla_sc.data(:,1))+10])


set(gca,'Tag','Sun_Scan');
grid on;
printfiles_report(gcf,Cal.dir_figs,'LineMode','Scaled');


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
                                           'date_range',datenum(2017,10,1):now,'process',0);  

events_cfg_op =getcfgs(tabla_dsp.data(:,1)',OP_config);    
events_cfg_chk=getcfgs(tabla_dsp.data(:,1)',ALT_config);    

% Command-Print 
 fprintf('\r\nDSP Analysis\r\n');  
 if size(tabla_dsp.data,1)<13      
    displaytable(tabla_dsp.data(:,[2:8 15:end])',tabla_dsp.events,12,'.4f',tabla_dsp.data_lbl([1:7 14:end]));
 else
    displaytable(tabla_dsp.data(:,[2:8 15:end]),tabla_dsp.data_lbl([2:7 14:end]),10,'.4f',tabla_dsp.events);
 end
 matrix2latex_ctable(tabla_dsp.data(:,15:end),fullfile(Cal.file_latex,['table_dsp_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                    'rowlabels',tabla_dsp.events,'columnlabels',tabla_dsp.data_lbl(14:end),...
                                    'alignment','c','resize',0.9);
                                
figure; set(gcf,'Tag',strcat('DSP_A1',Cal.brw_str{Cal.n_inst})); 
tabla_dsp.data(isnan(tabla_dsp.data(:,2)),:)=[];
l1=ploty(tabla_dsp.data(:,[1 15 17]),'s');
set(l1(1),'Marker','s','LineStyle','-','Color','k','MarkerFaceColor','k'); 
set(l1(2),'Marker','s','LineStyle','-','Color','g','MarkerFaceColor','g');
lim=get(gca,'XLim'); set(gca,'YLim',[0.3380 0.3430],'XLim',[lim(1)-10 lim(2)+10]); 
hold on; s(2)=stairs([events_cfg_chk.data(1,:) tabla_dsp.data(end,1)],[events_cfg_chk.data(7,:) events_cfg_chk.data(7,end)],'-','color','m','LineWidth',2);
         s(1)=stairs([events_cfg_op.data(1,:) tabla_dsp.data(end,1)],[events_cfg_op.data(7,:) events_cfg_op.data(7,end)],'-','color','b','LineWidth',2);
vl=vline_v(events_cfg_op.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3)); 
set(vl,'LineStyle','None'); set(findobj(gca,'Type','text'),'FontSize',7);    
title(sprintf('A1: %s',Cal.brw_name{Cal.n_inst}));
datetick('x',12,'KeepLimits','keepTicks'); ylabel('Ozone Absorption Coefficient (A1)'); grid
legendflex([l1;s(1);s(2)],{'A1 Quad','A1 Cubic','A1 ref (Op)','A1 ref (Chk)'}, ...
                           'anchor', {'sw','sw'}, 'buffer',[10 -5], 'nrow',2,'fontsize',8,'box','on');
set(findobj(gcf,'Type','text'),'BackgroundColor','w');
m183=nanmean(tabla_dsp.data(:,[15,17]));
hline(m183(1),'k',sprintf('quad=%.4f ',m183(1)))
hline(m183(2),'g',sprintf(' cubic=%.4f',m183(2)))

%%
figure
mq=mean_smooth_abs(tabla_dsp.data(:,1),tabla_dsp.data(:,[15]),90,0);
mc=mean_smooth_abs(tabla_dsp.data(:,1),tabla_dsp.data(:,[17]),90,0);
hold on
m183=nanmean(tabla_dsp.data(:,[15,17]));


plot(tabla_dsp.data(:,1),(mq(:,1)-m183(1))/0.001)
hold on
plot(tabla_dsp.data(:,1),(mc(:,1)-m183(2))/0.001)
legend({'cuadratic','cubic'})
plot(tabla_dsp.data(:,1),(tabla_dsp.data(:,15)-m183(1))/0.001,'k+')
hold on
plot(tabla_dsp.data(:,1),(tabla_dsp.data(:,17)-m183(2))/0.001,'bx')

datetick('x','mmm/yy','keepticks')
grid on
ylabel('Micromenter step')
title({'Brewer  #183 Ozone Abs coeff diff vs mean',sprintf('quad=%.4f cubic=%.4f',m183)})

set(gca,'YLim',[-2 2])
vline_v(event_info.dates,'k',event_info.labels);
hfill([-0.5,0.5],'grey')

%printfiles_report(gcf,Cal.dir_figs,'Width',14,'Height',7,'LineMode','Scaled');