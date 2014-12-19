% options_pub.outputDir=fullfile(pwd,'html');   options_pub.showCode=false; 
% close all; publish(fullfile(pwd,'Triad_comp.m'),options_pub);

%%  Brewer setup
clear all;
file_setup='triad_setup';
run(fullfile('.','cfg',file_setup));     % configuracion por defecto

Cal.Date.day0=datenum(2012,11,1);  
Cal.Date.dayend=datenum(2014,12,31); 
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend;
Cal.Date=Date;

Cal.file_save  = 'triad_comp_chk'; 
% Cal.file_latex = fullfile('..','latex'); 
% mkdir(Cal.file_latex);
% Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures'); 
% mkdir(Cal.dir_figs);

Cal.analyzed_brewer=1:3;
Cal.reference_brw=1:3;


%% ALL THE CONFIGURATIONS CONSTANTS HERE


osc_interval=[400,700,1000,1200];
Cal.brw=[157 183 185]; 

Cal.sl_c=[1,1,1];


% Generic
% Langley processing
airm_rang=[1.15 3.75]; 
% from condfiguration
%Fcorr={[0,0,0,0,0,0],[0,0,0,0,-8,0],[0,0,0,13,15,0]}; 
%
N_data=12; 
O3_std=2.5;
AOD_file=fullfile('.','130101_141231_Izana.lev15'); 
CLOUD_file=fullfile('..','cloudScreening.txt');

% Langley plots
ylim_brw={[1550 1650],[1500 1650],[1500 1650]}; 
ylim_dbs={[-50 50],[-50 50],[-50 50]};

grp_def=@(x) {year(x) weeknum(x)};




%% Update Brewer Summaries
load(Cal.file_save);
 for i=union(Cal.analyzed_brewer,Cal.reference_brw)
    
    processed_days=fix(cellfun(@(x) x(1,1),ozone_sum{i},'un',true));
    days_to_process=setdiff(Cal.Date.CALC_DAYS,processed_days);
    
    %ozone_raw{i}={};   hg{i}={};   ozone_sum{i}={};  ozone_raw0{i}={};  
    %config{i}={};      sl{i}={};   ozone_ds{i}={};   sl_cr{i}={};    

    Cal.Date.CALC_DAYS=days_to_process;
     
    [ozone,log_,missing_]=read_bdata(i,Cal);

    ozone_sum{i}=ozone.ozone_sum;
    config{i}=ozone.config;
    ozone_ds{i}=ozone.ozone_ds;
    ozone_raw{i}=ozone.raw;
    ozone_raw0{i}=ozone.raw0;
    sl{i}=ozone.sl;       % first calibration / bfiles
    sl_cr{i}=ozone.sl_cr; % recalc. with 2? configuration
 end
save(Cal.file_save);

%% Configs: Operative
for brw=1:3
    try
      events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,Cal.brw_config_files{brw,1});    
      fprintf('\nBrewer %s: Operative Config.\n',Cal.brw_name{brw});
      displaytable(cat(1,events_cfg_op.data(2:8,:),events_cfg_op.data(9,:)*10^8,events_cfg_op.data([10 17:end],:))',...
                         events_cfg_op.legend([2:10 17:end]),8,'.5g',cellstr(datestr(events_cfg_op.data(1,:),1))');
    catch exception
      fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{brw});
    end
end
%% Configs: Check
for brw=1:3
    try
      events_cfg_chk=getcfgs(Cal.Date.CALC_DAYS,Cal.brw_config_files{brw,2});    
      fprintf('\nBrewer %s: Second Config.\n',Cal.brw_name{brw});
      displaytable(cat(1,events_cfg_chk.data(2:8,:),events_cfg_chk.data(9,:)*10^8,events_cfg_chk.data([10 17:end],:))',...
                         events_cfg_chk.legend([2:10 17:end]),8,'.5g',cellstr(datestr(events_cfg_chk.data(1,:),1))');
    catch exception
      fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{brw});
    end
end

%% SL Report
close all;
for ii=union(Cal.analyzed_brewer,Cal.reference_brw)
    sl_mov_o{ii}={}; sl_median_o{ii}={}; sl_out_o{ii}={}; R6_o{ii}={};
    sl_mov_n{ii}={}; sl_median_n{ii}={}; sl_out_n{ii}={}; R6_n{ii}={};
% old instrumental constants
    [sl_mov_o{ii},sl_median_o{ii},sl_out_o{ii},R6_o{ii}]=sl_report_jday(ii,sl,Cal.brw_str,...
                               'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii},...
                               'hgflag',1,'fplot',0);
% new instrumental constants
    [sl_mov_n{ii},sl_median_n{ii},sl_out_n{ii},R6_n{ii}]=sl_report_jday(ii,sl_cr,Cal.brw_str,...
                               'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii},...
                               'hgflag',1,'fplot',0);
end
save(Cal.file_save,'-APPEND','sl_mov_o','sl_median_o','sl_out_o','R6_o',...
                             'sl_mov_n','sl_median_n','sl_out_n','R6_n');
                                                                          
%% Loading cfg. parameters
close all
[A,ETC,SL_B,SL_R,F_corr,SL_flag,cfg]=read_cal_config_new(config,Cal,{sl_median_o,sl_median_n});

%% Creating Summaries
for i=union(Cal.analyzed_brewer,Cal.reference_brw)
    cal{i}={}; summary_orig{i}={}; summary_orig_old{i}={};
    [cal{i},summary_orig{i},summary_orig_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,...
                                                                'flag_sl',1,'plot_sl',1);
    % filter correction
    [summary_old{i} summary{i}]=filter_corr(summary_orig,summary_orig_old,i,A,F_corr{i});
end
save(Cal.file_save,'-APPEND','A','ETC','F_corr','SL_B','SL_R','cfg','SL_flag',...
                             'summary_old','summary_orig_old','summary','summary_orig');

%% Comparacion Operativa.

[ref_op,ratio_ref_op]=join_summary(Cal,summary_old,Cal.reference_brw,Cal.analyzed_brewer,5);

% ratio_ref_plots
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref);
 
%% Comparacion Alternativa


[ref_alt,ratio_ref_alt]=join_summary(Cal,summary,Cal.reference_brw,Cal.analyzed_brewer,5);

% ratio_ref_plots
[f_hist_alt,f_ev_alt,f_sc_alt,f_smooth_alt]=ratio_ref_plots(Cal,ratio_ref_alt);


save(fullfile('.',Cal.file_save),'-APPEND','ref_op','ratio_ref_op',...
                                  'ref_alt','ratio_ref_alt')




%% Langley
close all

%% ---- langley from Indiv. Measurements ----
for brw=1:3
    [ozone_lgl{brw},cfg_indv,leg,ozone_lgl_sum{brw}] = langley_data_cell(ozone_raw{brw},ozone_ds{brw},config{brw});
end
   
%% ---- langley from Indiv. Measurements ---- 
for brw=1:3
    ozone_lgl_dep{brw}=langley_filter_lvl1(ozone_lgl{brw},'plots',0,...
                          'F_corr',F_corr{brw},'airmass',airm_rang,'O3_hday',O3_std,...
                          'AOD',AOD_file,'lgl_days',0,'plots',0);%...
%                        'Cloud',CLOUD_file);                   
    [brw_indv{brw} dbs_indv{brw} st_brw{brw} st_dbs{brw}] = langley_analys(ozone_lgl_dep,brw,Cal,...
                          'res_filt',1,'plot_flag',0);
end
save(fullfile('.',Cal.file_save),'-APPEND','ozone_lgl','cfg_indv','leg','ozone_lgl_sum',...
                                           'ozone_lgl_dep','brw_indv','st_brw','dbs_indv','st_dbs');
      
                                       
                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                       
                                       
                                       
                                       
                                       
%% R6 + ETC Langley
cfgs=1; 
for brw=1:3
    events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,Cal.brw_config_files{brw,1});    

    dates=cell2mat(events_raw{brw}(:,2)); indx=dates>=brw_indv{brw}(1,1) & dates<=brw_indv{brw}(end,1); 
    mixed_brw_indv=cat(1,brw_indv{brw}(:,[1 2],cfgs),brw_indv{brw}(:,[1 3],cfgs)); mixed_brw_indv=sortrows(mixed_brw_indv(~isnan(mixed_brw_indv(:,2)),:),1);
    mixed_dbs_indv=cat(1,dbs_indv{brw}(:,[1 2],cfgs),dbs_indv{brw}(:,[1 3],cfgs)); mixed_dbs_indv=sortrows(mixed_dbs_indv(~isnan(mixed_dbs_indv(:,2)),:),1);

    [m_brw,s_brw,n_brw,std_brw]=grpstats(mixed_brw_indv,grp_def(mixed_brw_indv(:,1)),{'mean','sem','numel','std'});
    [m_dbs,s_dbs,n_dbs,std_dbs]=grpstats(mixed_dbs_indv,grp_def(mixed_dbs_indv(:,1)),{'mean','sem','numel','std'});
    lmu=sortrows(m_brw,1); lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2)),1); lstd=sortrows(cat(2,m_brw(:,1),n_brw(:,2)),1);

    figure; set(gcf,'Tag',sprintf('Langley_R6_%s',Cal.brw_str{brw}));
    ha=tight_subplot(2,1,.08,.1,.075); hold all;
    axes(ha(1)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
    axes(ha(2)); set(gca,'box','on','YTickLabelMode','auto'); grid; hold on;

    axes(ha(1)); boundedline(gca,lmu(:,1),lmu(:,2),lsem(:,2),'-k',lmu(:,1),lmu(:,2),lstd(:,2),'-k',...
                                            'alpha' ,'transparency', 0.3);
    % Events
    if ~isempty(dates(indx))
       vl=vline(dates(indx),'r-'); set(vl,'LineWidth',1); 
    end
    % Config from icf
    idx=group_time(Cal.Date.CALC_DAYS',events_cfg_op.data(1,:));
    stairs(gca,cat(2,events_cfg_op.data(1,logical(unique(idx))),brw_indv{brw}(end,1,cfgs)),cat(2,events_cfg_op.data(8,unique(idx)),events_cfg_op.data(8,end)),'-g','LineWidth',2);
    ylabel('ETC (Brw method)');
    title(sprintf('Langley plot (%s - %s): %s, airmass range = [%.2f, %.2f]',...
              datestr(brw_indv{brw}(1,1,cfgs),28),datestr(brw_indv{brw}(end,1,cfgs),28),Cal.brw_name{brw},airm_rang)); 
                                          
    axes(ha(2)); 
    mmplotyy(SL_B.old(:,1),SL_B.old(:,brw+1),'gs',matadd(SL_R.old(:,brw+1),-SL_B.old(:,brw+1)),'r*');    
    ylabel('R6'); mmplotyy('R6 correction'); 
    % Events
    if ~isempty(dates(indx))
       vl=vline(dates(indx),'r-'); set(vl,'LineWidth',1); 
    end
    % Config from icf
    idx=group_time(Cal.Date.CALC_DAYS',events_cfg_op.data(1,:));
    s=findobj(gcf,'Type','axes');
    hold on; stairs(s(1),cat(2,events_cfg_op.data(1,logical(unique(idx))),brw_indv{brw}(end,1,cfgs)),...
                         cat(2,events_cfg_op.data(17,unique(idx)),events_cfg_op.data(17,end)),'-g','LineWidth',2);
    lg=legend('R6','R6 ref - R6 calc.','Location','NorthEast','Orientation','Horizontal'); set(lg,'FontSize',7);
    title(sprintf('SL R6 Ratio & SL Correction. Operative Cal. Constants.')); 
    datetick('x',6,'KeepTicks','KeepLimits');      
    set(ha(1),'XLim',[Cal.Date.CALC_DAYS(1)-10 Cal.Date.CALC_DAYS(end)+10]); linkprop(ha,'XLim');     
end 

%% Comparamos resultados (Diffs)
clc
for brw=1:3
    events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,Cal.brw_config_files{brw,1});    

    % difference am/pm con las dos confgs.
    fprintf('\nBrewer %s: AM vs PM\n',Cal.brw_name{brw});
    lgl_am_b(brw,:)=nanmedian(brw_indv{brw}(:,2,:)-brw_indv{brw}(:,3,:));
    lgl_am_d(brw,:)=nanmedian(dbs_indv{brw}(:,2,:)-dbs_indv{brw}(:,3,:));
    displaytable(cat(2,lgl_am_b(brw,:),lgl_am_d(brw,:)),{'Brw cfg1','Brw cfg2','Dbs cfg1','Dbs cfg2'},12,'.5g',{'AM vs PM'});

    % difference cfg1/cfg2 para am y pm
    fprintf('\nBrewer %s: CFG Op. vs CFG Chk\n',Cal.brw_name{brw});
    lgl_cfgs_b(brw,:)=nanmedian(brw_indv{brw}(:,2:3,1)-brw_indv{brw}(:,2:3,2));
    lgl_cfgs_d(brw,:)=nanmedian(dbs_indv{brw}(:,2:3,1)-dbs_indv{brw}(:,2:3,2));
    displaytable(cat(2,lgl_cfgs_b(brw,:),lgl_cfgs_d(brw,:)),{'Brw AM','Brw PM','Dbs AM','Dbs PM'},12,'.5g',{'CFG1 vs CFG2'});
    
    % tabla por eventos
    ampm_brw_op  =cat(1,brw_indv{brw}(:,[1 2],1),brw_indv{brw}(:,[1 3],1)); ampm_brw_op =sortrows(ampm_brw_op(~isnan(ampm_brw_op(:,2)),:),1);
    ampm_brw_chk =cat(1,brw_indv{brw}(:,[1 2],2),brw_indv{brw}(:,[1 3],2)); ampm_brw_chk=sortrows(ampm_brw_chk(~isnan(ampm_brw_chk(:,2)),:),1);
   
    ampm_dbs_op  =cat(1,dbs_indv{brw}(:,[1 2],1),dbs_indv{brw}(:,[1 3],1)); ampm_dbs_op =sortrows(ampm_dbs_op(~isnan(ampm_dbs_op(:,2)),:),1);
    ampm_dbs_chk =cat(1,dbs_indv{brw}(:,[1 2],2),dbs_indv{brw}(:,[1 3],2)); ampm_dbs_chk=sortrows(ampm_dbs_chk(~isnan(ampm_dbs_chk(:,2)),:),1);

    y=group_time(ampm_brw_op(:,1),events_cfg_op.data(1,:));

    fprintf('\nBrewer %s: CFG Op.\n',Cal.brw_name{brw});
    [m_b,s_b,n_b,std_b]=grpstats(ampm_brw_op,y,{'mean','sem','numel','std'});
    [m_d,s_d,n_d,std_d]=grpstats(ampm_dbs_op,y,{'mean','sem','numel','std'});
    lbls=Cal.events_raw{brw}(cell2mat(Cal.events_raw{brw}(:,2))>=events_cfg_op.data(1,1),3);
    displaytable(cat(2,m_b(:,2),s_b(:,2),n_b(:,2),m_d(:,2),s_d(:,2),n_d(:,2)),{'Brw','Brw sem','N','Dbs','Dbs sem','N'},...
        12,'.5g',lbls(unique(y))); 

    fprintf('\nBrewer %s: CFG Chk.\n',Cal.brw_name{brw});
    [m_b,s_b,n_b,std_b]=grpstats(ampm_brw_chk,y,{'mean','sem','numel','std'});
    [m_d,s_d,n_d,std_d]=grpstats(ampm_dbs_chk,y,{'mean','sem','numel','std'});
    displaytable(cat(2,m_b(:,2),s_b(:,2),n_b(:,2),m_d(:,2),s_d(:,2),n_d(:,2)),{'Brw','Brw sem','N','Dbs','Dbs sem','N'},...
        12,'.5g',lbls(unique(y))); 
end 
