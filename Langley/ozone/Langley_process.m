% options_pub=struct('outputDir',fullfile('..','..','html'),'showCode',false); 
% close all; publish(fullfile(pwd,'Langley157.m'),options_pub);

%% setup
clear all;
file_setup='calizo_setup';
run(fullfile('..','..',file_setup));     % configuracion por defecto
%eval('../Triad/load_data_157.m');
Cal.file_save  = fullfile('..','triad_2016_19')

try
    % to do not load Cal struct
    load(Cal.file_save,'A','ETC','F_corr','SL_B','SL_R','SL_corr_flag','cfg',...
                             'summary_old','summary_orig_old','summary','summary_orig','ozone_raw0','ozone_raw','ozone_ds','sl','sl_cr','config');
catch
    disp('new data');
end

Cal.Date.day0=datenum(2015,6,1);
Cal.Date.dayend=now;%datenum(2018,1,120);
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend; Cal.Date=Date;
%Cal.file_save  = fullfile('langley_2016_19')
%save(Cal.file_save,'A','ETC','F_corr','SL_B','SL_R','SL_corr_flag','cfg',...
%                             'summary_old','summary_orig_old','summary','summary_orig','ozone_raw0','ozone_raw','ozone_ds','sl','sl_cr','config')

%% Generic
for ii=1:3
    Cal.n_inst=ii
summ_old{Cal.n_inst}=summary_orig{Cal.n_inst};
summ{Cal.n_inst}=summary{Cal.n_inst};
% Langley processing
airm_rang={[1.15 4.0],[1.25 4.0],[1.15 4.0]};
Fcorr=F_corr;% {[0,0,0,0,0,0],[],[]}; Usaremos F's de la configuracion
N_data=12;
O3_std=2.5;
% generadas
AOD_file=fullfile('..','aod_2016_2019.lev15');
CLOUD_file=fullfile('..','cloudScreening.txt');

% Langley plots
ylim_brw={[1550 1650],[],[]};
ylim_dbs={[-50 50],[-50 50],[-50 50]};
grp_def=@(x) {year(x) weeknum(x)};


%% Configs: Operative
try
   events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,Cal.brw_config_files{Cal.n_inst,1});
   fprintf('\nBrewer %s: Operative Config.\n',Cal.brw_name{Cal.n_inst});
   displaytable(events_cfg_op.data(2:end,:),cellstr(datestr(events_cfg_op.data(1,:),1))',12,'.5g',events_cfg_op.legend(2:end));
catch exception
   fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{Cal.n_inst});
end

%% Configs: Check
try
   events_cfg_chk=getcfgs(Cal.Date.CALC_DAYS,Cal.brw_config_files{Cal.n_inst,2});
   fprintf('\nBrewer %s: Second Config.\n',Cal.brw_name{Cal.n_inst});
   displaytable(events_cfg_chk.data(2:end,:),cellstr(datestr(events_cfg_chk.data(1,:),1))',12,'.5g',events_cfg_chk.legend(2:end));
catch exception
   fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{Cal.n_inst});
end

                          

%% ---- langley from Indiv. Measurements ----
[ozone_lgl{Cal.n_inst},cfg_indv,leg,ozone_lgl_sum{Cal.n_inst}] = langley_data_cell(ozone_raw{Cal.n_inst},ozone_ds{Cal.n_inst},config{Cal.n_inst});
save(Cal.file_save,'-APPEND','ozone_lgl','cfg_indv','leg','ozone_lgl_sum');


%% ---- ND Filter :Langley
% Filters regression (aplicamos filtros de AOD y O3 stricto)
   cfgs=1;      % Operativa
   lgl_filt{Cal.n_inst}=langley_filter_lvl1(ozone_lgl{Cal.n_inst},'plots',0,...
                                    'airmass',[1.05 4.5],'O3_hday',2,'N_hday',N_data,...
                                 'AOD',AOD_file);%,'date_range',datenum(2014,1,[66 80]));

   brw_indv_{Cal.n_inst} = langley_analys_filter(lgl_filt,Cal.n_inst,...
                                         'res_filt',1,'plot_flag',0);

   nd0_=cat(1,brw_indv_{Cal.n_inst}(:,[1 2],cfgs),brw_indv_{Cal.n_inst}(:,[1 6],cfgs)); nd0=sortrows(nd0_,1);
   nd3_=cat(1,brw_indv_{Cal.n_inst}(:,[1 3],cfgs),brw_indv_{Cal.n_inst}(:,[1 7],cfgs)); nd3=sortrows(nd3_,1);
   nd4_=cat(1,brw_indv_{Cal.n_inst}(:,[1 4],cfgs),brw_indv_{Cal.n_inst}(:,[1 8],cfgs)); nd4=sortrows(nd4_,1);
   [m_brw,se_brw,s_brw,n_brw]=grpstats([nd0(:,1) nd3(:,2)-nd0(:,2),nd4(:,2)-nd0(:,2)],grp_def(nd0(:,1)),{'mean','sem','std','numel'});
   lmu{Cal.n_inst}=sortrows([m_brw,se_brw(:,2:end),s_brw(:,2:end),n_brw(:,2:end)]);
   %lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2:end)),1);
   brw_indv_filter{Cal.n_inst}=brw_indv_{Cal.n_inst};
  % save(Cal.file_save,'-APPEND','lgl_filt','brw_indv_filter');

   %% ---- langley from Indiv. Measurements ----

 [ozone_lgl_dep{Cal.n_inst},days_lgl{Cal.n_inst},days_all{Cal.n_inst}] =...
    langley_filter_lvl1(ozone_lgl{Cal.n_inst},'plots',0,...
    'F_corr',Fcorr{Cal.n_inst},'airmass',airm_rang{Cal.n_inst},'O3_hday',2.0,...
    'lgl_days',1,'plots',0,...
     'AOD',AOD_file,'Cloud',CLOUD_file);
                    %,'Cloud',CLOUD_file);

 [brw_indv{Cal.n_inst} dbs_indv{Cal.n_inst} st_brw{Cal.n_inst} st_dbs{Cal.n_inst}] = langley_analys(ozone_lgl_dep,Cal.n_inst,Cal,...
                       'res_filt',1,'plot_flag',0);
                   
 
                   
 %% salidas         
 t_days=array2table(days_lgl{Cal.n_inst}.data,'VariableNames',varname(days_lgl{Cal.n_inst}.labels),'Rownames',cellstr(datestr(days_lgl{Cal.n_inst}.data(:,1))));                  
 l_data_op=array2table([brw_indv{Cal.n_inst}(:,:,1),dbs_indv{Cal.n_inst}(:,2:end,1)],...
     'VariableNames',{'Diaj','F_AM','F_PM','S_AM','S_PM','D_AM','D_PM','DS_AM','DS_PM'},...
     'Rownames',cellstr(datestr(brw_indv{Cal.n_inst}(:,1,1))));  
 
 l_data_ch=array2table([brw_indv{Cal.n_inst}(:,:,1),dbs_indv{Cal.n_inst}(:,2:end,1)],'VariableNames',{'Diaj','F_AM','F_PM','S_AM','S_PM','D_AM','D_PM','DS_AM','DS_PM'},'Rownames',cellstr(datestr(brw_indv{Cal.n_inst}(:,1,2))));                  
 t_days_all=array2table(days_all{Cal.n_inst}.data,'VariableNames',varname(days_all{Cal.n_inst}.labels),'Rownames',cellstr(datestr(days_all{Cal.n_inst}.data(:,1))));                   

 %unimos langley y AOD/BSRN
 %lgl_table{Cal.n_inst}=join(t_days,l_data_op);   
 % all days table. 
 l_data_op.Properties.VariableNames{1}='Date'; 
 t_days_all.Properties.VariableNames{1}='Date'; 
 langley_table{Cal.n_inst}=outerjoin(l_data_op,t_days_all);   
 langley_table{Cal.n_inst}.Date=datetime(datestr(langley_table{Cal.n_inst}.Date_t_days_all));  
 langley_table{Cal.n_inst}=timetable2table(table2timetable(langley_table{Cal.n_inst}));
 writetable(langley_table{Cal.n_inst},'Langley_.xls','Sheet',Cal.brw_str{Cal.n_inst})

 write_langley(Cal.brw(Cal.n_inst),1619,brw_indv(Cal.n_inst),dbs_indv(Cal.n_inst));
 %save(Cal.file_save,'-APPEND','brw_indv','dbs_indv','st_brw','st_dbs','days_lgl','ozone_lgl_dep','langley_table');
 
end

% save
save(Cal.file_save,'-APPEND','lgl_filt','brw_indv_filter');
save(Cal.file_save,'-APPEND','brw_indv','dbs_indv','st_brw','st_dbs','days_lgl','ozone_lgl_dep','langley_table');
 

%% langley por eventos
% fprintf('Langley por eventos\r\n');
% event_info=getevents(Cal,'grp','events'); 
% cfgs=1;
% mixed_brw_indv=cat(1,brw_indv{Cal.n_inst}(:,[1 2],cfgs),brw_indv{Cal.n_inst}(:,[1 3],cfgs));
% mixed_brw_indv=sortrows(mixed_brw_indv(~isnan(mixed_brw_indv(:,2)),:),1);
% mixed_dbs_indv=cat(1,dbs_indv{Cal.n_inst}(:,[1 2],cfgs),dbs_indv{Cal.n_inst}(:,[1 3],cfgs));
% mixed_dbs_indv=sortrows(mixed_dbs_indv(~isnan(mixed_dbs_indv(:,2)),:),1);
% dates=cell2mat(events_raw{Cal.n_inst}(:,2)); indx=dates>=brw_indv{Cal.n_inst}(1,1) & dates<=brw_indv{Cal.n_inst}(end,1);
% 
% data_tab_brw=meanperiods(mixed_brw_indv, event_info); data_tab_dbs=meanperiods(mixed_dbs_indv, event_info);
% data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)*10)/10 round(data_tab_dbs.m(:,2)) round(data_tab_dbs.std(:,2)*10)/10 data_tab_dbs.N(:,2)];
% % Command-Print 
% t_events=display_table(cat(2,events_cfg_op.data(8,:)',data),{'ETC Op.','brw','std','dbs','std','N'},8,{'d','d','.1f','d','.1f','d'},data_tab_brw.evnts)
% 
% brw_lgl_op{Cal.n_inst}=mixed_brw_indv;
% dbs_lgl_op{Cal.n_inst}=mixed_dbs_indv;
% 
% 
% %save(Cal.file_save,'-APPEND','brw_lgl_op','dbs_lgl_op')
% 
% %% plot (Op. config)
%   cfgs=1; [a b]=fileparts(Cal.brw_config_files{Cal.n_inst,cfgs});
% 
%   figure; ha=tight_subplot(2,1,.08,.1,.075); hold all;
%   axes(ha(1)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
%   axes(ha(2)); set(gca,'box','on','YTickLabelMode','auto'); grid; hold on;
% 
%   plot(ha(1),brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,2,cfgs),'bd',brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,3,cfgs),'rd');
%   plot(ha(2),dbs_indv{Cal.n_inst}(:,1,cfgs),dbs_indv{Cal.n_inst}(:,2,cfgs),'bd',dbs_indv{Cal.n_inst}(:,1,cfgs),dbs_indv{Cal.n_inst}(:,3,cfgs),'rd');
%   set(findobj(ha,'Type','Line'),'MarkerSize',2);
%   legend({'AM','PM'},'Location','North','Orientation','Horizontal','FontSize',7);
% 
%  % Config from icf
%   idx=group_time(Cal.Date.CALC_DAYS',events_cfg_op.data(1,:));
%   events_cfg_op.data(1,logical(unique(idx)));
%   stairs(ha(1),cat(2,events_cfg_op.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs)),cat(2,events_cfg_op.data(8,unique(idx)),events_cfg_op.data(8,end)),'-k','LineWidth',2);
%   stairs(ha(2),cat(2,events_cfg_op.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs)),repmat(0*ones(1,length(unique(idx))+1),1),'-k','LineWidth',2);
% 
%   mixed_brw_indv=cat(1,brw_indv{Cal.n_inst}(:,[1 2],cfgs),brw_indv{Cal.n_inst}(:,[1 3],cfgs)); mixed_brw_indv=sortrows(mixed_brw_indv(~isnan(mixed_brw_indv(:,2)),:),1);
%   mixed_dbs_indv=cat(1,dbs_indv{Cal.n_inst}(:,[1 2],cfgs),dbs_indv{Cal.n_inst}(:,[1 3],cfgs)); mixed_dbs_indv=sortrows(mixed_dbs_indv(~isnan(mixed_dbs_indv(:,2)),:),1);
%   dates=cell2mat(events_raw{Cal.n_inst}(:,2)); indx=dates>=brw_indv{Cal.n_inst}(1,1) & dates<=brw_indv{Cal.n_inst}(end,1);
% 
% %  y=group_time(mixed_brw_indv(:,1),dates); % mean by events
%   [m_brw,s_brw,n_brw,std_brw]=grpstats(mixed_brw_indv,grp_def(mixed_brw_indv(:,1)),{'mean','sem','numel','std'});
%   [m_dbs,s_dbs,n_dbs,std_dbs]=grpstats(mixed_dbs_indv,grp_def(mixed_brw_indv(:,1)),{'mean','sem','numel','std'});
%   hold on; e1=errorbar(ha(1),m_brw(:,1),m_brw(:,2),s_brw(:,2),'d-.','MarkerFaceColor','g');
%            e2=errorbar(ha(2),m_dbs(:,1),m_dbs(:,2),s_dbs(:,2),'d-.','MarkerFaceColor','g');
%   set([e1 e2],'MarkerSize',6);  datetick('x',12,'keeplimits','keepticks');
%   title(ha(1),sprintf('Langley plot (%s - %s): %s, airmass range = [%.2f, %.2f], Config = %s',...
%               datestr(brw_indv{Cal.n_inst}(1,1,1),22),datestr(brw_indv{Cal.n_inst}(end,1,1),22),Cal.brw_name{Cal.n_inst},airm_rang,b));
%   ylabel(ha(1),'ETC (Brw method)','FontSize',8); ylabel(ha(2),'ETC corr. (Dbs method)','FontSize',8);
%   set(ha(1),'XLim',[Cal.Date.CALC_DAYS(1)-10 Cal.Date.CALC_DAYS(end)+10],...
%             'YLim',ylim_brw{Cal.n_inst}); linkprop(ha,'XLim'); set(ha(1),'YLim',ylim_brw{Cal.n_inst});
% 
% % Events
%   axes(ha(1)); vl=vline(dates(indx),'r-'); set(vl,'LineWidth',2);
%   axes(ha(2)); vl=vline(dates(indx),'r-'); set(vl,'LineWidth',2);
% 
%   snapnow
% 
% %%
% displaytable([m_brw(:,2),s_brw(:,2),m_dbs(:,2),s_dbs(:,2),n_brw(:,2)],...
%        {'ETC','std','ETC corr. (Dobson)','std (Dobson)','N'},13,'.0f',cellstr(datestr(m_brw(:,1),1)));
% tl_157=array2table([m_brw(:,2),s_brw(:,2),m_dbs(:,2),s_dbs(:,2),n_brw(:,2)],...
%     'VariableNames',str2name( {'ETC','std','ETC corr. (Dobson)','std (Dobson)','N'}),...
%     'RowNames',cellstr(datestr(m_brw(:,1),1)));
% tl_157.Date=datetime(datestr(m_brw(:,1),1))
% writetable(tl_157,'Langley_2016_2019.xls','Sheet',['Op_t_',Cal.brw_str{Cal.n_inst}],'WriteRownames',true);
% 
% %%
% fprintf('Langley por eventos\r\n');
% event_info=getevents(Cal,'grp','events'); 
% event_info.dates=events_cfg_op.data(1,:);
% event_info.labels=cellstr(datestr(events_cfg_op.data(1,:)));
% 
% data_tab_brw=meanperiods(mixed_brw_indv, event_info); 
% data_tab_dbs=meanperiods(mixed_dbs_indv, event_info);
% data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)*10)/10 round(data_tab_dbs.m(:,2)) round(data_tab_dbs.std(:,2)*10)/10 data_tab_dbs.N(:,2)];
% % Command-Print 
% displaytable(cat(2,events_cfg_op.data(8,:)',data),{'ETC Op.','brw','std','dbs','std','N'},8,{'d','d','.1f','d','.1f','d'},data_tab_brw.evnts);
% 
% %% Bounded plot
% figure; set(gcf,'Tag',sprintf('Lag%s_bound',Cal.brw_str{Cal.n_inst}));
% plot(brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,2,cfgs),'b.',...
%      brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,3,cfgs),'r.');
% datetick('x',12,'keeplimits','keepticks'); grid;
% [m_brw,s_brw,n_brw,std_brw]=grpstats(mixed_brw_indv,grp_def(mixed_brw_indv(:,1)),{'mean','sem','numel','std'});
% [m_dbs,s_dbs,n_dbs,std_dbs]=grpstats(mixed_dbs_indv,grp_def(mixed_dbs_indv(:,1)),{'mean','sem','numel','std'});
% lmu=sortrows(m_brw,1); lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2)),1); lstd=sortrows(cat(2,m_brw(:,1),n_brw(:,2)),1);
% boundedline(gca,lmu(:,1),lmu(:,2),lsem(:,2),'-k',lmu(:,1),lmu(:,2),lstd(:,2),'-k',...
%                                             'alpha' ,'transparency', 0.3);
% title(sprintf('Langley plot (%s - %s): %s, airmass range = [%.2f, %.2f]',...
%             datestr(brw_indv{Cal.n_inst}(1,1,cfgs),28),datestr(brw_indv{Cal.n_inst}(end,1,cfgs),28),Cal.brw_name{Cal.n_inst},airm_rang));
% ylabel('ETC (Brw method)'); set(gca,'XLim',[brw_indv{Cal.n_inst}(1,1,cfgs)-5 brw_indv{Cal.n_inst}(end,1,cfgs)+5]);
% % Events
% vl=vline(dates(indx),'r-'); set(vl,'LineWidth',1);
% % Config from icf
%   idx=group_time(Cal.Date.CALC_DAYS',events_cfg_op.data(1,:));
%   stairs(gca,cat(2,events_cfg_op.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs)),cat(2,events_cfg_op.data(8,unique(idx)),events_cfg_op.data(8,end)),'-g','LineWidth',2);
% 
% %% R6 + ETC Langley
% %   figure; set(gcf,'Tag',sprintf('Langley_R6_%s',Cal.brw_str{Cal.n_inst}));
% %   ha=tight_subplot(2,1,.08,.1,.075); hold all;
% %   axes(ha(1)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
% %   axes(ha(2)); set(gca,'box','on','YTickLabelMode','auto'); grid; hold on;
% % 
% %   axes(ha(1)); boundedline(gca,lmu(:,1),lmu(:,2),lsem(:,2),'-k',lmu(:,1),lmu(:,2),lstd(:,2),'-k',...
% %                                             'alpha' ,'transparency', 0.3);
% %   % Events
% %   vl=vline(dates(indx),'r-'); set(vl,'LineWidth',1);
% %   % Config from icf
% %   idx=group_time(Cal.Date.CALC_DAYS',events_cfg_op.data(1,:));
% %   stairs(gca,cat(2,events_cfg_op.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs)),cat(2,events_cfg_op.data(8,unique(idx)),events_cfg_op.data(8,end)),'-g','LineWidth',2);
% %   ylabel('ETC (Brw method)');
% %   title(sprintf('Langley plot (%s - %s): %s, airmass range = [%.2f, %.2f]',...
% %             datestr(brw_indv{Cal.n_inst}(1,1,cfgs),28),datestr(brw_indv{Cal.n_inst}(end,1,cfgs),28),Cal.brw_name{Cal.n_inst},airm_rang));
% % 
% %   axes(ha(2));
% %   mmplotyy(SL_B.old(:,1),SL_B.old(:,Cal.n_inst+1),'gs',matadd(SL_R.old(:,Cal.n_inst+1),-SL_B.old(:,Cal.n_inst+1)),'r*');
% %   ylabel('R6'); mmplotyy('R6 correction');
% %   % Events
% %   vl=vline(dates(indx),'r-'); set(vl,'LineWidth',1);
% %   % Config from icf
% %   idx=group_time(Cal.Date.CALC_DAYS',events_cfg_op.data(1,:));
% %   s=findobj(gcf,'Type','axes');
% %   hold on; stairs(s(1),cat(2,events_cfg_op.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs)),...
% %                        cat(2,events_cfg_op.data(17,unique(idx)),events_cfg_op.data(17,end)),'-g','LineWidth',2);
% %   lg=legend('R6','R6 ref - R6 calc.','Location','NorthEast','Orientation','Horizontal'); set(lg,'FontSize',7);
% %   title(sprintf('SL R6 Ratio & SL Correction. Operative Cal. Constants.')); 
% %   datetick('x',6,'KeepTicks','KeepLimits');      
% %   set(ha(1),'XLim',[Cal.Date.CALC_DAYS(1)-10 Cal.Date.CALC_DAYS(end)+10]); linkprop(ha,'XLim');     
% 
% %% plot  (Chk. config)
%   cfgs=2; [a b]=fileparts(Cal.brw_config_files{Cal.n_inst,cfgs});
% 
%   figure; ha=tight_subplot(2,1,.08,.1,.075); hold all;
%   axes(ha(1)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
%   axes(ha(2)); set(gca,'box','on','YTickLabelMode','auto'); grid; hold on;
% 
%   plot(ha(1),brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,2,cfgs),'bd',brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,3,cfgs),'rd'); 
%   plot(ha(2),dbs_indv{Cal.n_inst}(:,1,cfgs),dbs_indv{Cal.n_inst}(:,2,cfgs),'bd',dbs_indv{Cal.n_inst}(:,1,cfgs),dbs_indv{Cal.n_inst}(:,3,cfgs),'rd'); 
%   set(findobj(ha,'Type','Line'),'MarkerSize',2);
%   lg=legend({'AM','PM'},'Location','North','Orientation','Horizontal','FontSize',7);   
%  
%  % Config from icf
%   idx=group_time(Cal.Date.CALC_DAYS',events_cfg_chk.data(1,:));
%   %events_cfg_chk.data(1,logical(unique(idx)))
%   %xx=cat(2,events_cfg_chk.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs))
%   %yy=cat(2,events_cfg_chk.data(8,unique(idx)),events_cfg_chk.data(8,end))
%   stairs(ha(1),cat(2,events_cfg_chk.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs)),'-k','LineWidth',2);
%   stairs(ha(2),cat(2,events_cfg_chk.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs)),repmat(0*ones(1,length(unique(idx))+1),1),'-k','LineWidth',2);
%  
%   mixed_brw_indv=cat(1,brw_indv{Cal.n_inst}(:,[1 2],cfgs),brw_indv{Cal.n_inst}(:,[1 3],cfgs));
%   mixed_brw_indv=sortrows(mixed_brw_indv(~isnan(mixed_brw_indv(:,2)),:),1);
%   mixed_dbs_indv=cat(1,dbs_indv{Cal.n_inst}(:,[1 2],cfgs),dbs_indv{Cal.n_inst}(:,[1 3],cfgs)); 
%   mixed_dbs_indv=sortrows(mixed_dbs_indv(~isnan(mixed_dbs_indv(:,2)),:),1);
%   dates=cell2mat(events_raw{Cal.n_inst}(:,2)); indx=dates>=brw_indv{Cal.n_inst}(1,1) & dates<=brw_indv{Cal.n_inst}(end,1); 
% 
%   brw_lgl_ck{Cal.n_inst}=mixed_brw_indv;
%   dbs_lgl_ck{Cal.n_inst}=mixed_dbs_indv;
% 
% 
%   %save(Cal.file_save,'-APPEND','brw_lgl_ck','dbs_lgl_ck')
% 
%   
%   
%   
% %  y=group_time(mixed_brw_indv(:,1),dates); % mean by events
%   [m_brw,s_brw,n_brw,std_brw]=grpstats(mixed_brw_indv,grp_def(mixed_brw_indv(:,1)),{'mean','sem','numel','std'});
%   [m_dbs,s_dbs,n_dbs,std_dbs]=grpstats(mixed_dbs_indv,grp_def(mixed_dbs_indv(:,1)),{'mean','sem','numel','std'});
%   hold on; e1=errorbar(ha(1),m_brw(:,1),m_brw(:,2),s_brw(:,2),'d-.','MarkerFaceColor','g');
%            set(e1,'MarkerSize',6);
%            e2=errorbar(ha(2),m_dbs(:,1),m_dbs(:,2),s_dbs(:,2),'d-.','MarkerFaceColor','g');
%            set(e2,'MarkerSize',6);         
%   datetick('x',12,'keeplimits','keepticks');
%   title(ha(1),sprintf('Langley plot (%s - %s): %s, airmass range = [%.2f, %.2f], Config = %s',...
%               datestr(brw_indv{Cal.n_inst}(1,1,1),22),datestr(brw_indv{Cal.n_inst}(end,1,1),22),Cal.brw_name{Cal.n_inst},airm_rang,b)); 
%   ylabel(ha(1),'ETC (Brw method)','FontSize',8);
%   ylabel(ha(2),'ETC corr. (Dbs method)','FontSize',8);
%   set(ha(1),'XLim',[Cal.Date.CALC_DAYS(1)-10 Cal.Date.CALC_DAYS(end)+10],...
%             'YLim',ylim_brw{Cal.n_inst}); linkprop(ha,'XLim');
%         set(ha(1),'YLim',ylim_brw{Cal.n_inst});
%          set(ha(2),'YLim',[-35 35]);
% 
%   % Events
%   axes(ha(1)); vl=vline(dates(indx),'r-'); set(vl,'LineWidth',2); 
%   axes(ha(2)); vl=vline(dates(indx),'r-'); set(vl,'LineWidth',2); 
%    
%   snapnow
% %%
% fprintf('Langley por eventos\r\n');
% event_info=getevents(Cal,'grp','events'); 
% event_info.dates=events_cfg_op.data(1,:);
% event_info.labels=cellstr(datestr(events_cfg_op.data(1,:)));
% 
% 
% data_tab_brw=meanperiods(mixed_brw_indv, event_info); 
% data_tab_dbs=meanperiods(mixed_dbs_indv, event_info);
% data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)*10)/10 round(data_tab_dbs.m(:,2)) round(data_tab_dbs.std(:,2)*10)/10 data_tab_dbs.N(:,2)];
% % Command-Print 
% % try
% % displaytable(cat(2,events_cfg_chk.data(8,:)',data)',{'ETC Op.','brw','std','dbs','std','N'},8,{'d','d','.1f','d','.1f','d'},data_tab_brw.evnts);
% % catch
% %     disp('error');
% % end
% 
% display_table(cat(2,events_cfg_chk.data(8,:)',data),{'ETC Op.','brw','std','dbs','std','N'},8,{'d','d','.1f','d','.1f','d'},data_tab_brw.evnts);
% 
% 
% %%
% displaytable([m_brw(:,2),s_brw(:,2),m_dbs(:,2),s_dbs(:,2),n_brw(:,2)],...
%        {'ETC','std','ETC corr. (Dobson)','std (Dobson)','N'},13,'.0f',cellstr(datestr(m_brw(:,1),1)));
% tl_157=array2table([m_brw(:,2),s_brw(:,2),m_dbs(:,2),s_dbs(:,2),n_brw(:,2)],...
%     'VariableNames',str2name( {'ETC','std','ETC corr. (Dobson)','std (Dobson)','N'}),...
%     'RowNames',cellstr(datestr(m_brw(:,1),1)));
% tl_157.Date=datetime(datestr(m_brw(:,1),1))
% writetable(tl_157,'Langley_2017.xls','Sheet',['Al_t_',Cal.brw_str{Cal.n_inst}],'WriteRownames',true);
% %% Bounded plot
% figure; set(gcf,'Tag',sprintf('Lag%s_bound',Cal.brw_str{Cal.n_inst}));
% plot(brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,2,cfgs),'b.',...
%      brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,3,cfgs),'r.'); 
% datetick('x',12,'keeplimits','keepticks'); grid;
% [m_brw,s_brw,n_brw,std_brw]=grpstats(mixed_brw_indv,grp_def(mixed_brw_indv(:,1)),{'mean','sem','numel','std'});
% [m_dbs,s_dbs,n_dbs,std_dbs]=grpstats(mixed_dbs_indv,grp_def(mixed_dbs_indv(:,1)),{'mean','sem','numel','std'});
% lmu=sortrows(m_brw,1); lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2)),1); lstd=sortrows(cat(2,m_brw(:,1),n_brw(:,2)),1);
% boundedline(gca,lmu(:,1),lmu(:,2),lsem(:,2),'-k',lmu(:,1),lmu(:,2),lstd(:,2),'-k',...
%                                             'alpha' ,'transparency', 0.3);
% title(sprintf('Langley plot (%s - %s): %s, airmass range = [%.2f, %.2f]',...
%             datestr(brw_indv{Cal.n_inst}(1,1,cfgs),28),datestr(brw_indv{Cal.n_inst}(end,1,cfgs),28),Cal.brw_name{Cal.n_inst},airm_rang)); 
% ylabel('ETC (Brw method)'); set(gca,'XLim',[brw_indv{Cal.n_inst}(1,1,cfgs)-5 brw_indv{Cal.n_inst}(end,1,cfgs)+5]);
% % Events
% vl=vline(dates(indx),'r-'); set(vl,'LineWidth',1); 
% % Config from icf
% idx=group_time(Cal.Date.CALC_DAYS',events_cfg_chk.data(1,:));
% stairs(gca,cat(2,events_cfg_chk.data(1,logical(unique(idx))),brw_indv{Cal.n_inst}(end,1,cfgs)),cat(2,events_cfg_chk.data(8,unique(idx)),events_cfg_chk.data(8,end)),'-g','LineWidth',2);
% 
% %% Comparamos las configuraciones posibles: 
% 
% %% ploteo de los residuos
% % Op. Constants
% brw_op=cell2mat(cellfun(@(x) x(~isnan(x(:,4,1)),4,1),st_brw{Cal.n_inst}.r,'UniformOutput',0));
% dbs_op=cell2mat(cellfun(@(x) x(~isnan(x(:,4,1)),4,1),st_dbs{Cal.n_inst}.r,'UniformOutput',0));
% 
% % Chk. Constants
% brw_chk=cell2mat(cellfun(@(x) x(~isnan(x(:,4,2)),4,2),st_brw{Cal.n_inst}.r,'UniformOutput',0));
% dbs_chk=cell2mat(cellfun(@(x) x(~isnan(x(:,4,2)),4,2),st_dbs{Cal.n_inst}.r,'UniformOutput',0));
% 
% figure; set(gcf,'Tag',sprintf('Lag%s_residuals',Cal.brw_str{Cal.n_inst}));
% aux={brw_op,dbs_op,brw_chk,dbs_chk};
% tex=nhist(aux,'samebins','box','smooth','ylabel','Pdf (Langley fit residuals)');
% legend({'Brw Op. (DT=26)','Dbs. Op. (DT=23)','Brw Chk. (DT=26)','Dbs. Chk. (DT=23)'},'location','NorthEast','FontSize',7);
% title(Cal.brw_name{Cal.n_inst}); set(findobj(gca,'Type','Line'),'Marker','None','LineWidth',1)
% 
% %% Rsquare
% figure; set(gcf,'Tag',sprintf('Lag%s_Rsquare',Cal.brw_str{Cal.n_inst}));
% plot(st_brw{Cal.n_inst}.rs(:,1,1),st_brw{Cal.n_inst}.rs(:,2,1),'.m',...
%      st_brw{Cal.n_inst}.rs(:,1,1),st_brw{Cal.n_inst}.rs(:,2,2),'.g');
% ylabel('Rsquare'); title(sprintf('Brewer %s',Cal.brw_name{Cal.n_inst}));
% datetick('x',19,'keeplimits','keepticks'); set(gca,'YLim',[0.9995 1]); grid;
% legend({'BrewerMth: DT 26 ns','BrewerMth: DT 23 ns'},'Location','SouthEast','Orientation','Horizontal');
%   
% %% Comparamos resultados (Diffs)
% clc
%  % difference am/pm con las dos confgs.
%  fprintf('\nBrewer %s: AM vs PM\n',Cal.brw_name{Cal.n_inst});
%  lgl_am_b(Cal.n_inst,:)=nanmedian(brw_indv{Cal.n_inst}(:,2,:)-brw_indv{Cal.n_inst}(:,3,:));
%  lgl_am_d(Cal.n_inst,:)=nanmedian(dbs_indv{Cal.n_inst}(:,2,:)-dbs_indv{Cal.n_inst}(:,3,:));
%  displaytable(cat(2,lgl_am_b(Cal.n_inst,:),lgl_am_d(Cal.n_inst,:)),{'Brw cfg1','Brw cfg2','Dbs cfg1','Dbs cfg2'},12,'.5g',{'AM vs PM'});
% 
%  % difference cfg1/cfg2 para am y pm
%  fprintf('\nBrewer %s: CFG Op. vs CFG Chk\n',Cal.brw_name{Cal.n_inst});
%  lgl_cfgs_b(Cal.n_inst,:)=nanmedian(brw_indv{Cal.n_inst}(:,2:3,1)-brw_indv{Cal.n_inst}(:,2:3,2));
%  lgl_cfgs_d(Cal.n_inst,:)=nanmedian(dbs_indv{Cal.n_inst}(:,2:3,1)-dbs_indv{Cal.n_inst}(:,2:3,2));
%  displaytable(cat(2,lgl_cfgs_b(Cal.n_inst,:),lgl_cfgs_d(Cal.n_inst,:)),{'Brw AM','Brw PM','Dbs AM','Dbs PM'},12,'.5g',{'CFG1 vs CFG2'});
%  
%  % tabla por eventos
%  ampm_brw_op  =cat(1,brw_indv{Cal.n_inst}(:,[1 2],1),brw_indv{Cal.n_inst}(:,[1 3],1)); ampm_brw_op =sortrows(ampm_brw_op(~isnan(ampm_brw_op(:,2)),:),1);
%  ampm_brw_chk =cat(1,brw_indv{Cal.n_inst}(:,[1 2],2),brw_indv{Cal.n_inst}(:,[1 3],2)); ampm_brw_chk=sortrows(ampm_brw_chk(~isnan(ampm_brw_chk(:,2)),:),1);
% 
%  ampm_dbs_op  =cat(1,dbs_indv{Cal.n_inst}(:,[1 2],1),dbs_indv{Cal.n_inst}(:,[1 3],1)); ampm_dbs_op =sortrows(ampm_dbs_op(~isnan(ampm_dbs_op(:,2)),:),1);
%  ampm_dbs_chk =cat(1,dbs_indv{Cal.n_inst}(:,[1 2],2),dbs_indv{Cal.n_inst}(:,[1 3],2)); ampm_dbs_chk=sortrows(ampm_dbs_chk(~isnan(ampm_dbs_chk(:,2)),:),1);
% 
%  y=group_time(ampm_brw_op(:,1),events_cfg_op.data(1,:));
% 
%  fprintf('\nBrewer %s: CFG Op.\n',Cal.brw_name{Cal.n_inst});
%  [m_b,s_b,n_b,std_b]=grpstats(ampm_brw_op,y,{'mean','sem','numel','std'});
%  [m_d,s_d,n_d,std_d]=grpstats(ampm_dbs_op,y,{'mean','sem','numel','std'});
%  lbls=Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op.data(1,1),3);
%  displaytable(cat(2,m_b(:,2),s_b(:,2),n_b(:,2),m_d(:,2),s_d(:,2),n_d(:,2)),{'Brw','Brw sem','N','Dbs','Dbs sem','N'},...
%      12,'.5g',lbls(unique(y))); 
% 
%  fprintf('\nBrewer %s: CFG Chk.\n',Cal.brw_name{Cal.n_inst});
%  [m_b,s_b,n_b,std_b]=grpstats(ampm_brw_chk,y,{'mean','sem','numel','std'});
%  [m_d,s_d,n_d,std_d]=grpstats(ampm_dbs_chk,y,{'mean','sem','numel','std'});
%  displaytable(cat(2,m_b(:,2),s_b(:,2),n_b(:,2),m_d(:,2),s_d(:,2),n_d(:,2)),{'Brw','Brw sem','N','Dbs','Dbs sem','N'},...
%      12,'.5g',lbls(unique(y))); 
