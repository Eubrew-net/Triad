% options_pub=struct('outputDir',fullfile('..','html'),'showCode',true)
% close all; publish(fullfile(pwd,'load_data.m'),options_pub);

%%  Brewer setup
clear all
run(fullfile('..','read_config_'))

%Cal.Date.day0=datenum(2018,2,1);
%Cal.Date.dayend=datenum(2018,3,10);
Cal.Date.day0=datenum(2019,7,10);
Cal.Date.dayend=datenum(2019,7,15);
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend; Cal.Date=Date;

brewer=[157,183,185];

% READ Brewer Summaries
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

 for ii=Cal.n_ref
    
    sl_mov_o{ii}={}; sl_median_o{ii}={}; sl_out_o{ii}={}; R6_o{ii}={};
    sl_mov_n{ii}={}; sl_median_n{ii}={}; sl_out_n{ii}={}; R6_n{ii}={};
    % old instrumental constants
    [sl_mov_o{ii},sl_median_o{ii},sl_out_o{ii},R6_o{ii}]=sl_report_jday(ii,sl,Cal.brw_str,...
                               'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii}(:,1:5),...
                               'hgflag',1,'fplot',0);
 end
 
 [A,ETC,SL_B,SL_R,F_corr,SL_corr_flag,cfg]=read_cal_config_new(config,Cal,{sl_median_o,sl_median_n});

% Data recalculation for summaries  and individual observations
for i=Cal.n_ref
    cal{i}={}; summary_orig{i}={}; summary_orig_old{i}={};
    [cal{i},summary_orig{i},summary_orig_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,...
                                          'flag_sl',1,'plot_sl',0,'flag_sl_corr',SL_corr_flag);
    % filter correction
    [summary_old{i} summary{i}]=filter_corr(summary_orig,summary_orig_old,i,A,F_corr{i});
end


 b1=1;
 b2=3;
 b3=2;
 
 clear rat
 nds=size(summary_old{b1},1);
 for i=1:nds
    rat(i,1)=summary_old{b1}(i,1);
    rat(i,2)=summary_old{b1}(i,4);
    rat(i,3)=summary_old{b1}(i,8);
    rat(i,15)=summary_old{b1}(i,3);
    rat(i,4)=summary{b1}(i,8);
    % primer ratio
    [m,ri]=min(abs(summary_old{b2}(:,1)-summary_old{b1}(i,1)));
    if abs(summary_old{b2}(ri,1)-summary_old{b1}(i,1))<3./24/60; %la diff de tiempo entre las medidas < 1 min
        rat(i,5)=summary_old{b2}(ri,8);
        rat(i,6)=summary{b2}(ri,8);
        rat(i,9)=rat(i,3)/rat(i,5);
        rat(i,10)=rat(i,4)/rat(i,6);
    else
        rat(i,5)=NaN;
        rat(i,6)=NaN;
        rat(i,9)=NaN;
        rat(i,10)=NaN;
    end
    % segund ratio
    [m,ri]=min(abs(summary_old{b3}(:,1)-summary_old{b1}(i,1)));
    if abs(summary_old{b3}(ri,1)-summary_old{b1}(i,1))<3./24/60; %la diff de tiempo entre las medidas < 1 min
        rat(i,7)=summary_old{b3}(ri,8);
        rat(i,8)=summary{b3}(ri,8);
        rat(i,11)=rat(i,3)/rat(i,7);
        rat(i,12)=rat(i,4)/rat(i,8);
    else
        rat(i,7)=NaN;
        rat(i,8)=NaN;
        rat(i,11)=NaN;
        rat(i,12)=NaN;
    end
    
    if ~isnan(rat(i,5)) & ~isnan(rat(i,7))
        rat(i,13)=rat(i,5)/rat(i,7);
    else
        rat(i,13)=NaN;
    end   

    if ~isnan(rat(i,6)) & ~isnan(rat(i,8))
        rat(i,14)=rat(i,6)/rat(i,8);
    else
        rat(i,14)=NaN;
    end   
    
 end
 
figure
plot(rat(:,1),rat(:,9),'.')
%ylim([prctile([rat(:,7);rat(:,11)],1),prctile([rat(:,7);rat(:,11)],99)])
hold
plot(rat(:,1),rat(:,11),'.')
plot(rat(:,1),rat(:,13),'.')
hline(nanmedian(rat(:,9)))
hline(nanmedian(rat(:,11)))
hline(nanmedian(rat(:,13)))
legend(sprintf("ms9(%d)/ms9(%d)",brewer(b1),brewer(b2)),sprintf("ms9(%d)/ms9(%d)",brewer(b1),brewer(b3)),sprintf("ms9(%d)/ms9(%d)",brewer(b2),brewer(b3)))
% set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',1,'keeplimits','keepticks')
title(sprintf("ms9 ratios %d, %d and %d Op. (%s,%s)",brewer(b1),brewer(b2),brewer(b3),datestr(Date.CALC_DAYS(1)),datestr(Date.CALC_DAYS(end))));
xlabel("Date")
ylabel("ms9 ratio between brewers")

figure
plot(rat(:,1),rat(:,10),'.')
%ylim([prctile([rat(:,8);rat(:,12)],1),prctile([rat(:,8);rat(:,12)],99)])
hold
plot(rat(:,1),rat(:,12),'.')
plot(rat(:,1),rat(:,14),'.')
hline(nanmedian(rat(:,10)))
hline(nanmedian(rat(:,12)))
hline(nanmedian(rat(:,14)))
legend(sprintf("ms9(%d)/ms9(%d)",brewer(b1),brewer(b2)),sprintf("ms9(%d)/ms9(%d)",brewer(b1),brewer(b3)),sprintf("ms9(%d)/ms9(%d)",brewer(b2),brewer(b3)))
datetick('x',1,'keeplimits','keepticks')
title(sprintf("ms9 ratios %d, %d and %d Alt. (%s,%s)",brewer(b1),brewer(b2),brewer(b3),datestr(Date.CALC_DAYS(1)),datestr(Date.CALC_DAYS(end))));
xlabel("Date")
ylabel("ms9 ratio between brewers")

figure
plot(rat(:,2),rat(:,13),'.')
hold
xlabel("Temperature")
ylabel(sprintf("ms9(%d)/ms9(%d)",brewer(b2),brewer(b3)))
[fitt stat]=robustfit(rat(:,2),rat(:,13));
xl = xlim;
yr = fitt(1) + fitt(2) * xl;
plot(xl,yr)
yl = ylim;

figure
plot(rat(:,15),rat(:,13),'.')
hold
xlabel("Air mass")
ylabel(sprintf("ms9(%d)/ms9(%d)",brewer(b2),brewer(b3)))
[fitm stat]=robustfit(rat(:,15),rat(:,13));
xl = xlim;
yr = fitm(1) + fitm(2) * xl;
plot(xl,yr)
yl = ylim;

figure
plot(rat(:,1),rat(:,13),'.')
hold
plot(rat(:,1),rat(:,13)-fitt(2) * rat(:,15),'.')
plot(rat(:,1),rat(:,13)-fitm(2) * rat(:,15),'.')
xlabel("Date")
ylabel("ms9 ratio between brewers")
legend(sprintf("ms9(%d)/ms9(%d) no corregido",brewer(b2),brewer(b3)),sprintf("ms9(%d)/ms9(%d) corregido temp",brewer(b2),brewer(b3)),sprintf("ms9(%d)/ms9(%d) corregido m.o.",brewer(b2),brewer(b3)))




