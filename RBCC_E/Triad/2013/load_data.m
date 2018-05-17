% options_pub.outputDir=fullfile(pwd,'..','html');   options_pub.showCode=false; 
% publish(fullfile(pwd,'load_data.m'),options_pub);

%%  Brewer setup
clear all;
file_setup='calizo2013_setup';
run(fullfile('..','..',file_setup));     % configuracion por defecto

Cal.Date.day0=datenum(2013,1,1);  Cal.Date.dayend=datenum(2013,12,31); 
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend; Cal.Date=Date;

Cal.file_save  = 'triad_comp2013'; 
Cal.file_latex = fullfile('..','latex'); mkdir(Cal.file_latex);
Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures'); mkdir(Cal.dir_figs);

% % Op. config
% Cal.brw_config_files{1,1}=fullfile(Cal.path_root,'campaigns','Triad_maintenance','cfg','config157.cfg');
% Cal.brw_config_files{2,1}=fullfile(Cal.path_root,'campaigns','Triad_maintenance','cfg','config183.cfg');
% Cal.brw_config_files{3,1}=fullfile(Cal.path_root,'campaigns','Triad_maintenance','cfg','config185.cfg');
% % Chk. config
% Cal.brw_config_files{1,2}=fullfile(Cal.path_root,'campaigns','Triad_maintenance','cfg','config157_a.cfg');
% Cal.brw_config_files{2,2}=fullfile(Cal.path_root,'campaigns','Triad_maintenance','cfg','config183_a.cfg');
% Cal.brw_config_files{3,2}=fullfile(Cal.path_root,'campaigns','Triad_maintenance','cfg','config185_a.cfg');
 
%% READ Brewer Summaries
Cal.n_ref=Cal.n_ref;
 for i=Cal.n_ref
    ozone_raw{i}={};   hg{i}={};   ozone_sum{i}={};  ozone_raw0{i}={};  
    config{i}={};      sl{i}={};   ozone_ds{i}={};   sl_cr{i}={};    

    [ozone,log_,missing_]=read_bdata(i,Cal);

    % depuramos datos (ver incidencias en config. matrix)
    ozone=dep_data(Cal.incidences_text{i},ozone);
    
    ozone_sum{i}=ozone.ozone_sum;
    config{i}=ozone.config;
    ozone_ds{i}=ozone.ozone_ds;
    ozone_raw{i}=ozone.raw;
    ozone_raw0{i}=ozone.raw0;
    sl{i}=ozone.sl;       % first calibration / bfiles
    sl_cr{i}=ozone.sl_cr; % recalc. with 2? configuration
 end
save(Cal.file_save);

%% Configs
for i=Cal.n_ref
    %% Operative
    OP_config=Cal.brw_config_files{i,1};
    try
      [a b c]=fileparts(OP_config);
      if ~strcmpi(strcat(b,c),sprintf('config%d.cfg',Cal.brw(i)))
         fprintf(strcat(1,'\rCUIDADO!! Puede que las configuraciones cargadas no sean las esperadas\n',...
                          '(Expected: %s, Loading: %s)\n'),...
                 sprintf('config%d.cfg',Cal.brw(i)),strcat(b,c)); 
      end
      fprintf('\nBrewer %s: Operative Config.\n',Cal.brw_name{i});
      events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,OP_config);    
      displaytable(events_cfg_op.data(2:end,:),cellstr(datestr(events_cfg_op.data(1,:),1))',12,'.5g',events_cfg_op.legend(2:end));
    catch exception
      fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{i});
    end
    %% Check
    ALT_config=Cal.brw_config_files{i,2};
    try
       [a b c]=fileparts(ALT_config);
       if ~strcmpi(strcat(b,c),sprintf('config%d_a.cfg',Cal.brw(i)))
          fprintf(strcat(1,'\rCUIDADO!! Puede que las configuraciones cargadas no sean las esperadas\n',...
                           '(Expected: %s, Loading: %s)\n'),...
                  sprintf('config%d_a.cfg',Cal.brw(i)),strcat(b,c)); 
       end    
       events_cfg_chk=getcfgs(Cal.Date.CALC_DAYS,ALT_config);    
       fprintf('\nBrewer %s: Second Config.\n',Cal.brw_name{i});
       displaytable(events_cfg_chk.data(2:end,:),cellstr(datestr(events_cfg_chk.data(1,:),1))',12,'.5g',events_cfg_chk.legend(2:end));
    catch exception
       fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{i});
    end
end

%% SL Report
close all;
for ii=Cal.n_ref
    sl_mov_o{ii}={}; sl_median_o{ii}={}; sl_out_o{ii}={}; R6_o{ii}={};
    sl_mov_n{ii}={}; sl_median_n{ii}={}; sl_out_n{ii}={}; R6_n{ii}={};
% old instrumental constants
    [sl_mov_o{ii},sl_median_o{ii},sl_out_o{ii},R6_o{ii}]=sl_report_jday(ii,sl,Cal.brw_str,...
                               'outlier_flag',0,'diaj_flag',0,'events_raw',events_raw{ii},...
                               'hgflag',1,'fplot',1);
 % Imprimimos valores por eventos  
    fprintf('SL means, Op. config (by events). Brewer %s\r\n',Cal.brw_name{ii}); Cal.n_inst=ii;
    event_info=getevents(Cal,'grp','events'); data_tab=meanperiods(sl_median_o{ii}(:,[1 2 4]), event_info);
    displaytable(cat(2,data_tab.m(:,2),data_tab.std(:,2),data_tab.m(:,3),data_tab.std(:,3),data_tab.N(:,1)),...
                 {'R6','std','R5','std','N'},15,'.4f',data_tab.evnts);

% new instrumental constants
    [sl_mov_n{ii},sl_median_n{ii},sl_out_n{ii},R6_n{ii}]=sl_report_jday(ii,sl_cr,Cal.brw_str,...
                               'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii},...
                               'hgflag',1,'fplot',0);
% Imprimimos valores por eventos  
    fprintf('SL means, Chk. config (by events). Brewer %s\r\n',Cal.brw_name{ii}); Cal.n_inst=ii;
    event_info=getevents(Cal,'grp','events'); data_tab=meanperiods(sl_median_n{ii}(:,[1 2 4]), event_info);
    displaytable(cat(2,data_tab.m(:,2),data_tab.std(:,2),data_tab.m(:,3),data_tab.std(:,3),data_tab.N(:,1)),...
                 {'R6','std','R5','std','N'},15,'.4f',data_tab.evnts);

end
save(Cal.file_save,'-APPEND','sl_mov_o','sl_median_o','sl_out_o','R6_o',...
                             'sl_mov_n','sl_median_n','sl_out_n','R6_n');
                                                                          
%% SL report plots
ix=sort(findobj('Tag','SL_R6_report')); Width=15; Height=7; str=Cal.brw_str(Cal.n_ref);
for i=1:length(ix)
    figure(ix(i)); 
    datetick('x',19,'keeplimits','keepticks'); th=rotateticklabel(gca,45); set(th,'FontSize',9); 
    set(findobj(gcf,'Type','text'),'FontSize',8); set(gca,'XTicklabel','','FontSize',10); xlabel(''); 
    set(findobj(gca,'Marker','.'),'MarkerSize',5); set(findobj(gca,'Marker','o','-or','Marker','s'),'MarkerSize',4)
    if i~=2 
       set(findobj(gcf,'Tag','legend'),'FontSize',8,'Location','SouthWest','LineWidth',0.3);
    else
       set(findobj(gcf,'Tag','legend'),'FontSize',8,'Location','NorthWest','LineWidth',0.3);
    end
    % Config from Operative icf
    OP_config=Cal.brw_config_files{i,1};
    events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,OP_config);    
    idx=group_time(R6_n{i}(:,1),events_cfg_op.data(1,:));
    stairs(gca,R6_n{i}(logical(idx),1),events_cfg_op.data(17,idx),'-b','LineWidth',3);

    printfiles_report(gcf,Cal.dir_figs,'aux_pattern',str(i),'Width',Width,'Height',Height,'LineMode','Scaled');
end

%% Creating Summaries
close all
[A,ETC,SL_B,SL_R,F_corr,SL_corr_flag,cfg]=read_cal_config_new(config,Cal,{sl_median_o,sl_median_n});
% SL_corr_flag.old(~isnan(SL_corr_flag.old(:,3)),3)=1;

% Data recalculation for summaries  and individual observations
for i=Cal.n_ref
    cal{i}={}; summary_orig{i}={}; summary_orig_old{i}={};
    [cal{i},summary_orig{i},summary_orig_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,...
                                                             'flag_sl',1,'flag_sl_corr',SL_corr_flag,'plot_sl',1);
    % filter correction
    [summary_old{i} summary{i}]=filter_corr(summary_orig,summary_orig_old,i,A,F_corr{i});
end
write_summary(Cal.brw(1:3),Cal.Date.cal_year,summary_old,summary,SL_R,SL_B);
save(Cal.file_save,'-APPEND','A','ETC','F_corr','SL_B','SL_R','cfg',...
                             'summary_old','summary_orig_old','summary','summary_orig');
ref=[];
for ii=Cal.n_ref
    if Cal.sl_c(ii)
       med=summary_old{ii}(:,[1 12]); meds=summary_old{ii}(:,[1 13]);
    else
       med=summary_old{ii}(:,[1 6]); meds=summary_old{ii}(:,[1 7]);
    end
    % redondeamos la medida cada 5 minutos
    TSYNC=5;
    time=([fix(med(:,1)*24*60/TSYNC)/24/60*TSYNC,med(:,1)]);  
    med(:,1)= time(:,1); 
    ref=scan_join(ref,med);
end 
figure; set(gcf,'Tag','ozone');
ploty(ref,'.'); grid
legend(gca,Cal.brw_name{Cal.n_ref},'Location','NorthEast','Orientation','Horizontal'); 
datetick('x',6,'Keeplimits','Keepticks'); title('DS Ozone'); ylabel('Ozone (DU)');

%% Comparacion Operativa.
% close all

reference_brw=[1 2 3]; analyzed_brewer=[1 2 3];  
osc_interval=[400,700,1000,1200];
Cal.analyzed_brewer=analyzed_brewer; Cal.sl_c=[0,0,0];

Cal.brw=[157 183 185];
[ref,ratio_ref]=join_summary(Cal,summary_old,reference_brw,analyzed_brewer,10);

% osc table
[ratio_osc_table,osc_matrix,osc_stats]=osc_table(Cal,ratio_ref,osc_interval);
set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid

% ratio_ref_plots
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref);
 
%% Comparacion Alternativa
% close all

reference_brw=[1 2 3]; analyzed_brewer=[1 2 3];  
osc_interval=[400,700,1000,1200];
Cal.analyzed_brewer=analyzed_brewer; Cal.sl_c=[0,0,0];

[ref,ratio_ref]=join_summary(Cal,summary,reference_brw,analyzed_brewer,10);

% osc table
[ratio_osc_table,osc_matrix,osc_stats]=osc_table(Cal,ratio_ref,osc_interval);
set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid

% ratio_ref_plots
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref);
