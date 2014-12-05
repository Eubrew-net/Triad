% options_pub.outputDir=fullfile(pwd,'html');   options_pub.showCode=false; 
% close all; publish(fullfile(pwd,'Triad_comp.m'),options_pub);

%%  Brewer setup
clear all;
file_setup='triad_setup';
run(fullfile('.','cfg',file_setup));     % configuracion por defecto

Cal.Date.day0=datenum(2013,12,15);  
Cal.Date.dayend=datenum(2014,1,15); 
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend;
Cal.Date=Date;

Cal.file_save  = 'triad_comp_'; 
% Cal.file_latex = fullfile('..','latex'); 
% mkdir(Cal.file_latex);
% Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures'); 
% mkdir(Cal.dir_figs);

Cal.analyzed_brewer=1:3;
Cal.reference_brw=1:3;
 
%% READ Brewer Summaries
 for i=union(Cal.analyzed_brewer,Cal.reference_brw)
    ozone_raw{i}={};   hg{i}={};   ozone_sum{i}={};  ozone_raw0{i}={};  
    config{i}={};      sl{i}={};   ozone_ds{i}={};   sl_cr{i}={};    

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
[A,ETC,SL_B,SL_R,F_corr,cfg]=read_cal_config_new(config,Cal,{sl_median_o,sl_median_n});

%% Creating Summaries
for i=union(Cal.analyzed_brewer,Cal.reference_brw)
    cal{i}={}; summary_orig{i}={}; summary_orig_old{i}={};
    [cal{i},summary_orig{i},summary_orig_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,...
                                                                'flag_sl',Cal.sl_c(i),'plot_sl',1);
    % filter correction
    [summary_old{i} summary{i}]=filter_corr(summary_orig,summary_orig_old,i,A,F_corr{i});
end
save(Cal.file_save,'-APPEND','A','ETC','F_corr','SL_B','SL_R','cfg',...
                             'summary_old','summary_orig_old','summary','summary_orig');

%% Comparacion Operativa.
close all

osc_interval=[400,700,1000,1200];

Cal.brw=[157 183 185];
[ref,ratio_ref]=join_summary(Cal,summary_old,Cal.reference_brw,Cal.analyzed_brewer,5);

% ratio_ref_plots
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref);
 
%% Comparacion Alternativa
close all

osc_interval=[400,700,1000,1200];

[ref,ratio_ref]=join_summary(Cal,summary,Cal.reference_brw,Cal.analyzed_brewer,5);

% ratio_ref_plots
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref);
