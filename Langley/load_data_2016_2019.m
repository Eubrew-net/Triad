% options_pub.outputDir=fullfile(pwd,'report');   options_pub.showCode=true;
% publish(fullfile(pwd,'load_data.m'),options_pub);

%%  Brewer setup
clear all;
file_setup='calizo_setup';
run(fullfile('..',file_setup));     % configuracion por defecto

Cal.Date.day0=datenum(2015,6,1);
Cal.Date.dayend=now; % to K&Z calibration
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend; Cal.Date=Date;

Cal.file_save  = 'triad_2016_19';
Cal.file_latex = fullfile('..','latex'); 
Cal.dir_tables = fullfile('..','latex','tables')
mkdir(Cal.file_latex);
Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures');
mkdir(Cal.dir_figs);
mkdir(Cal.dir_tables);
%load(Cal.file_save);
if exist(Cal.file_save,'file')
     save(Cal.file_save,'-append','Cal','events_raw'); % sobreescribe la configuracion
     load(Cal.file_save);
else
    disp('clean');
end



% READ Brewer Summaries
Cal.n_ref=Cal.n_ref;
 for i=1:Cal.n_brw
    ozone_raw{i}={};   hg{i}={};   ozone_sum{i}={};  ozone_raw0{i}={};
    config{i}={};      sl{i}={};   ozone_ds{i}={};   sl_cr{i}={};

    if i<=3
        [ozone,log_,missing_]=read_bdata(i,Cal);
    else
        [ozone,log_,missing_]=read_bdata(i,Cal,'/Users/aredondas/CODE/iberonesia/RBCC_E/2016/bdata');
    end
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
save(Cal.file_save,'ozone_sum','config','ozone_ds','ozone_raw','ozone_raw0','sl','sl_cr','Cal','hg');
%%
%load(Cal.file_save,'ozone_sum','config','ozone_ds','ozone_raw','ozone_raw0','sl','sl_cr')
%%
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
      op_cfg{i}=display_table(events_cfg_op.data(2:end,:),cellstr(datestr(events_cfg_op.data(1,:),1))',12,'.5g',events_cfg_op.legend(2:end));
    catch exception
      fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{i});
    end
    matrix2latex_config(events_cfg_op.data(2:end,:),fullfile(Cal.dir_tables,['Op_config_',Cal.brw_str{i},'.tex']),...
                    'rowlabels',events_cfg_op.legend(2:end),'columnlabels',cellstr(datestr(events_cfg_op.data(1,:),1))',...
                    'size','footnotesize');
    
    
    
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
       alt_cfg{i}=display_table(events_cfg_chk.data(2:end,:),cellstr(datestr(events_cfg_chk.data(1,:),1))',12,'.5g',events_cfg_chk.legend(2:end));
       
           matrix2latex_config(events_cfg_chk.data(2:end,:),fullfile(Cal.dir_tables,['Op_config_',Cal.brw_str{i},'.tex']),...
                    'rowlabels',events_cfg_chk.legend(2:end),'columnlabels',cellstr(datestr(events_cfg_chk.data(1,:),1))',...
                    'size','footnotesize');
       
       
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
                               'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii},...
                               'hgflag',1,'fplot',1);
    % Imprimimos valores por eventos
    fprintf('SL means, Op. config (by events). Brewer %s\r\n',Cal.brw_name{ii}); Cal.n_inst=ii;
    event_info=getevents(Cal,'grp','events'); 
    data_tab=meanperiods(sl_median_o{ii}(:,[1 2 4]), event_info);
    events_SL{ii}=display_table(cat(2,round(data_tab.m(:,2)),round(data_tab.std(:,2),1),round(data_tab.m(:,3)),round(data_tab.std(:,3),1),data_tab.N(:,1)),...
                 {'R6','r6_std','R5','r5_std','N'},15,'.4f',data_tab.evnts);
    events_SL{ii}.Date=datetime(datevec(fix(event_info.dates)));
    events_SL{ii}.Label=event_info.labels';
    movevars(events_SL{ii},'Date','Before','R6');
    events_SL{ii}.R6_ref= op_cfg{ii}{16,:}';
    movevars(events_SL{ii},'R6_ref','Before','R6')
    
    
    writetable(events_SL{ii},'R6_events.xls','Sheet',Cal.brw_str{ii})
    % SL report plots
    
  
     ix=maxf(findobj('Tag','SL_R6_report'));
     figure(ix);
    
     Width=15; Height=7;
     str=Cal.brw_str(Cal.n_ref);
     axis tight
     datetick('x',12,'keeplimits','keepticks');
     th=rotateticklabel(gca,45); set(th,'FontSize',9);
     set(findobj(gcf,'Type','text'),'FontSize',8); 
     set(gca,'XTicklabel','','FontSize',10); xlabel('');
     set(findobj(gca,'Marker','.'),'MarkerSize',5);
     set(findobj(gca,'Marker','o','-or','Marker','s'),'MarkerSize',4)
     set(findobj(gcf,'Tag','legend'),'FontSize',8,'Location','NorthWest','LineWidth',0.3);
    
     % Config from Operative icf
     OP_config=Cal.brw_config_files{ii,1};
     events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,OP_config);
     idx=group_time(R6_o{ii}(:,1),events_cfg_op.data(1,:));
     stairs(gca,R6_o{ii}(logical(idx),1),events_cfg_op.data(17,idx),'-b','LineWidth',3);
     set(ix,'Tag',['SL_R6_report_',Cal.brw_str{ii}])
     printfiles_report(gcf,Cal.dir_figs,'aux_pattern' ,str(ii),'Width',Width*3,'Height',Height*1.5,'LineMode','Scaled');
     snapnow; 
   
           
    disp(ii)
end
save(Cal.file_save,'-APPEND','sl_mov_o','sl_median_o','sl_out_o','R6_o',...
                             'sl_mov_n','sl_median_n','sl_out_n','R6_n');
                         
%% Creating Summaries
close all
[A,ETC,SL_B,SL_R,F_corr,SL_corr_flag,cfg]=read_cal_config_new(config,Cal,{sl_median_o,sl_median_n});

% Data recalculation for summaries  and individual observations
for i=Cal.n_ref
    cal{i}={}; summary_orig{i}={}; summary_orig_old{i}={};
    [cal{i},summary_orig{i},summary_orig_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,...
                                          'flag_sl',1,'plot_sl',1,'flag_sl_corr',SL_corr_flag);
    % filter correction
    [summary_old{i} summary{i}]=filter_corr(summary_orig,summary_orig_old,i,A,F_corr{i});
end
t_sum=write_summary_cfg((1:3),20162019,summary_old,summary,SL_R,SL_B,A,ETC);
save(Cal.file_save,'-APPEND','A','ETC','F_corr','SL_B','SL_R','SL_corr_flag','cfg',...
                             'summary_old','summary_orig_old','summary','summary_orig','t_sum');
%% Outliers? Los detectamos
ref=[];
for ii=Cal.n_ref
    med=summary_old{ii}(:,[1 6]); meds=summary_old{ii}(:,[1 7]);
    TSYNC=5;   time=([fix(med(:,1)*24*60/TSYNC)/24/60*TSYNC,med(:,1)]);
    med(:,1)= time(:,1);    ref=scan_join(ref,med);
end
figure; set(gcf,'Tag','ozone');
ploty(ref,'.'); grid
legend(gca,Cal.brw_name{Cal.n_ref},'Location','NorthEast','Orientation','Horizontal');
datetick('x',6,'Keeplimits','Keepticks'); title('DS Ozone'); ylabel('Ozone (DU)');

% %% Comparacion Operativa.
% % close all
% 
% reference_brw=1:3; analyzed_brewer=[1 2 3];
% osc_interval=[400,700,1000,1200];
% Cal.analyzed_brewer=analyzed_brewer;
% Cal.sl_c=[0,0,0];
% 
% Cal.brw=[157 183 185];
% % before, during and after huelva
% 
% % BEFORE Iza?a
% [ref,ratio_ref]=join_summary(Cal,summary_old,reference_brw,analyzed_brewer,5,'date_range',datenum(2016,1,[1,175]));
% % ratio_ref_plots
% %outliers
% Cal.analyzed_brewer=[1,2,3];
% Cal.reference_brewer=[1,2,3];
% for jj=1:length(Cal.analyzed_brewer)
% [ratio_ref(:,jj+1),b,c]=deleteoutliers(ratio_ref(:,jj+1),0.05,1);
%  disp(datestr(ratio_ref(b,1)))
%  disp(c)
% end
% [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',0);
% 
% [ratio_osc_table,osc_matrix]=osc_table(Cal,ratio_ref,osc_interval);
% set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% 
% 
% %%  DURING Iza?a
% sum_x=summary_old;
% %sum_x{2}=summary{2};
% [ref,ratio_ref]=join_summary(Cal,sum_x,1:2,1:2,5,'date_range',datenum(2016,1,[175,210]));
% % ratio_ref_plots
% Cal.analyzed_brewer=[1,2];
% Cal.reference_brewer=[1,2];
% 
% %outliers
% for jj=1:length(Cal.analyzed_brewer)
% [ratio_ref(:,jj+1),b,c]=deleteoutliers(ratio_ref(:,jj+1),0.05,1);
%  disp(datestr(ratio_ref(b,1)))
%  disp(c)
% end
% 
% [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
% 
% [ratio_osc_table,osc_matrix]=osc_table(Cal,ratio_ref,osc_interval);
% set(gcf,'Tag','osc_box_plot');
% set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% %%
% %%  DURING2
% sum_x=summary;
% [ref,ratio_ref]=join_summary(Cal,sum_x,[1:3],[1:3],5,'date_range',datenum(2016,1,[220,240]));
% % ratio_ref_plots
% Cal.analyzed_brewer=[1,2,3];
% Cal.reference_brewer=[1];
% 
% %outliers
% for jj=1:length(Cal.analyzed_brewer)
% [ratio_ref(:,jj+1),b,c]=deleteoutliers(ratio_ref(:,jj+1),0.05,1);
%  disp(datestr(ratio_ref(b,1)))
%  disp(c)
% end
% 
% [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
% 
% [ratio_osc_table,osc_matrix]=osc_table(Cal,ratio_ref,osc_interval);
% set(gcf,'Tag','osc_box_plot');
% set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% %% AFTER Iza?a
% [ref,ratio_ref]=join_summary(Cal,summary_old,1:3,1:3,5,'date_range',datenum(2016,1,[210,366]));
% % ratio_ref_plots
% Cal.analyzed_brewer=[1,2,3];
% Cal.reference_brewer=[1,2,3];
% %outliers
% for jj=1:length(Cal.analyzed_brewer)
% [ratio_ref(:,jj+1),b,c]=deleteoutliers(ratio_ref(:,jj+1),0.05,1);
%  disp(datestr(ratio_ref(b,1)))
%  disp(c)
% end
% 
% [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',0);
% 
% [ratio_osc_table,osc_matrix]=osc_table(Cal,ratio_ref,osc_interval);
% set(gcf,'Tag','osc_box_plot');
% set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% %% ATMOZ Iza?a
% [ref,ratio_ref]=join_summary(Cal,summary_old,1:3,1:3,5,'date_range',datenum(2016,1,[240,280]));
% % ratio_ref_plots
% Cal.analyzed_brewer=[1,2,3];
% Cal.reference_brewer=[1,2,3];
% %outliers
% for jj=1:length(Cal.analyzed_brewer)
% [ratio_ref(:,jj+1),b,c]=deleteoutliers(ratio_ref(:,jj+1),0.05,1);
%  disp(datestr(ratio_ref(b,1)))
%  disp(c)
% end
% 
% [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',0);
% 
% [ratio_osc_table,osc_matrix]=osc_table(Cal,ratio_ref,osc_interval);
% set(gcf,'Tag','osc_box_plot');
% set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% 
% 
% %% Comparacion Alternativa
% close all
% 
% reference_brw=[1 2 3]; analyzed_brewer=[1 2 3];
% osc_interval=[400,700,1000,1200];
% Cal.analyzed_brewer=analyzed_brewer;
% Cal.sl_c=[0,0,0];
% 
% 
% % BEFORE IZa?a
% [ref,ratio_ref_alt]=join_summary(Cal,summary,reference_brw,analyzed_brewer,5,'date_range',datenum(2016,1,[1,175]));
% % ratio_ref_plots
% for jj=1:length(Cal.analyzed_brewer)
% [ratio_ref_alt(:,jj+1),b,c]=deleteoutliers(ratio_ref_alt(:,jj+1),0.05,1);
%  disp(datestr(ratio_ref_alt(b,1)))
%  disp(c)
% end
% 
% [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref_alt,'plot_smooth',0);
% [ratio_osc_table,osc_matrix]=osc_table(Cal,ratio_ref_alt,osc_interval);
% set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% 
% %%  DURING Iza?a
% [ref,ratio_ref_al]=join_summary(Cal,summary,1:3,1:3,5,'date_range',datenum(2016,1,[250,273]));
% % ratio_ref_plots
% Cal.analyzed_brewer=[1,2,3];
% Cal.reference_brewer=[1,2,3];
% [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref_al,'plot_smooth',0);
% [ratio_osc_table,osc_matrix]=osc_table(Cal,ratio_ref_al,osc_interval);
%  set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% 
% 
% %% AFTER HUELVA
% [ref,ratio_ref_aln]=join_summary(Cal,summary,1:3,1:3,5,'date_range',datenum(2016,1,[210,366]));
% % ratio_ref_plots
% Cal.analyzed_brewer=[1,2,3];
% Cal.reference_brewer=[1,2,3];
% [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref_aln,'plot_smooth',0);
% 
% 
% [ratio_osc_table,osc_matrix]=osc_table(Cal,ratio_ref_aln,osc_interval);
%  
%  
%  
%  set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% 
% % [ref,ratio_ref]=join_summary(Cal,summary,reference_brw,analyzed_brewer,10);
% % % osc table
% % [ratio_osc_table,osc_matrix,osc_stats]=osc_table(Cal,ratio_ref,osc_interval);
% % set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% % % ratio_ref_plots
% % [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',0);
