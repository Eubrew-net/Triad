% options_pub.outputDir=fullfile(pwd,'html'); options_pub.showCode=false;
% close all; publish(fullfile(pwd,'proccesing.m'),options_pub);

%% Processing things: Brewer setup
clear all;
file_setup='triad_setup';
run(fullfile('.','cfg',file_setup));     % configuracion por defecto

Cal.file_latex = fullfile('.','latex'); mkdir(Cal.file_latex);
Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures'); mkdir(Cal.dir_figs);

Cal.brw=[157 183 185];
Cal.analyzed_brewer=1:3; Cal.reference_brw=1:3;

Cal.reference_brw=1;
%% Ozone deviations
year_period=2005:2014
Cal.sl_c=[1,1,1];
Cal.path_root='C:\CODE\iberonesia\RBCC_E\Triad'
% loading data (summaries)
summ_op=cell(1,length(Cal.brw));
summ_alt=cell(1,length(Cal.brw));
for brw=1:length(Cal.brw)    
    for yr=1:length(year_period) 
        b=year_period(yr)
        %dir_summ=fullfile(Cal.path_root,'..',num2str(year_period(yr)),'Triad','Langley');
        dir_summ=fullfile(Cal.path_root,num2str(year_period(yr)))
        file_op=dir(fullfile(dir_summ,strcat('summary_old_Brw',num2str(Cal.brw(brw)),'*')));
        a=load(fullfile(dir_summ,file_op.name));
         a1=import(fullfile(dir_summ,file_op.name));
        %summ_op{brw}=[summ_op{brw},a]
        %summ_op{brw}=[summ_op{brw};load(fullfile(dir_summ,file_op.name))]
        
        summ_op{brw}=cat(1,summ_op{brw},load(fullfile(dir_summ,file_op.name)));

        %file=dir(fullfile(dir_summ,strcat('summary_Brw',num2str(Cal.brw(brw)),'*')));
        %summ_alt{brw}=cat(1,summ_alt{brw},load(fullfile(dir_summ,file.name)));
    end    
end

% ploteo: Comparacion Operativa.
[ref_op,ratio_ref_op]=join_summary(Cal,summ_op,Cal.reference_brw,Cal.analyzed_brewer,5,...
                                                 'date_range',[datenum(2009,1,1) datenum(2015,12,31)]);
%                                                  'date_range',[datenum(2014,6,15) datenum(2014,8,17)]);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref_op,'plot_smooth',0);

% Otro
[a b]=grpstats(ratio_ref_op,{year(ratio_ref_op(:,1)),weeknum(ratio_ref_op(:,1))},{@(x) nanmean(x,1),@(x) nanstd(x,0,1)});
figure;  errorbar(a(:,1),a(:,2),b(:,2),'sk','MarkerFaceColor','k');
hold on; errorbar(a(:,1),a(:,3),b(:,3),'sb','MarkerFaceColor','b');
         errorbar(a(:,1),a(:,4),b(:,4),'sr','MarkerFaceColor','r');
set(gca,'Ylim',[-1.5 1.5]); datetick('x',12,'KeepLimits','keepTicks'); 
title('Ozone Deviation to the Triad Mean'); ylabel('Ozone Deviation [%]'); 
set(findobj(gcf,'Marker','s'),'MarkerSize',3); grid         
legend(cellfun(@(x) strcat('Brw#',x),mmcellstr(sprintf('%03d|',Cal.brw)),'UniformOutput',0),...
                    'Location','North','Orientation','Horizontal');

figure; hold all
plot(ref_op(:,1),ref_op(:,2),'k.',...
     ref_op(:,1),ref_op(:,3),'g.',...
     ref_op(:,1),ref_op(:,4),'r.');
set(findobj(gcf,'Marker','.'),'MarkerSize',5); grid
lg=legend('IZO#157','IZO#183','IZO#185','Location','South','Orientation','Horizontal'); 
set(lg,'FontSize',8,'FontWeight','Demi');
datetick('x',12,'Keeplimits','Keepticks'); title('Izaña, Canaries'); ylabel('Ozone Total Content [DU]');
 
%%
[a b]=grpstats(ratio_ref_op,{year(ratio_ref_op(:,1)),weeknum(ratio_ref_op(:,1))},{@(x) nanmean(x,1),@(x) nanstd(x,0,1)});
figure; set(gcf,'Tag','Triad_deviation'); ha=tight_subplot(2,1,.08,.1,.075); 

axes(ha(1)); 
plot(ref_op(:,1),ref_op(:,2),'k.',...
     ref_op(:,1),ref_op(:,3),'r.',...
     ref_op(:,1),ref_op(:,4),'g.');
set(gca,'YLim',[200 400]);
ylabel('Total Ozone [DU]'); title('RBCC-E Brewer Triad, Izana, Canary Islands'); 
set(gca,'FontSize',8); set(findobj(gcf,'Marker','.'),'MarkerSize',3); grid; box on;
datetick('x',12,'keeplimits','keepticks');
 
axes(ha(2)); errorbar(a(:,1),a(:,2),b(:,2),'sk','MarkerFaceColor','k');
hold on;     errorbar(a(:,1),a(:,3),b(:,3),'sr','MarkerFaceColor','r');
             errorbar(a(:,1),a(:,4),b(:,4),'sg','MarkerFaceColor','g');
set(gca,'Ylim',[-1.5 1.5]); datetick('x',12,'KeepLimits','keepTicks'); 
title('Ozone Deviation to the Triad Mean'); ylabel('Ozone Deviation [%]'); 
set(findobj(gcf,'Marker','s'),'MarkerSize',3); grid         
lg=legend(cellfun(@(x) strcat('Brw#',x),mmcellstr(sprintf('%03d|',Cal.brw)),'UniformOutput',0),...
                    'Location','SouthEast','Orientation','Horizontal'); set(lg,'Fontsize',6);

xL=get(ha(1),'XTick');
set(ha(1),'XTickLabel',[],'XLim',[xL(1)-10,xL(end)+10]);    linkprop(ha,'XLim');
                
%%
for brw=1:length(Cal.brw)    
    figure; set(gcf,'Tag',sprintf('R6_%s',Cal.brw_str{brw}));
    mmplotyy(summ_op{brw}(:,1),summ_op{brw}(:,end),'sg',...
             matadd(summ_op{brw}(:,end-1),-summ_op{brw}(:,end)),'r.'); hold on;
    set(findobj(gca,'Marker','s'),'MarkerSize',4);
    st=stairs(summ_op{brw}(:,1),summ_op{brw}(:,end-1),'-'); set(st,'LineWidth',2);
    title(sprintf('%s: R6 and SL correction (red). Izaña, Canaries',Cal.brw_name{brw}));
    ylabel('MS9'); a=mmplotyy('SL corr = SL ref - SL calc.');  set(a,'Color','r');
    datetick('x',12,'Keeplimits','Keepticks'); grid
end

%% Langley results
year_period=2009:2015;

% loading data (Langley results)
lang=cell(1,length(Cal.brw));
for brw=1:length(Cal.brw)    
    for yr=1:length(year_period)        
        dir_summ=fullfile(Cal.path_root,'..',num2str(year_period(yr)),'Triad','Langley','ozone');
        
        file_op=dir(fullfile(dir_summ,strcat('Langley_Brw',num2str(Cal.brw(brw)),'*')));
        lang{brw}=cat(1,lang{brw},load(fullfile(dir_summ,file_op.name)));
    end
end

% ploteo: Langley
ylim_brw={[1525 1675],[1500 1780],[1350 1650]}; 
grp_def=@(x) {year(x) weeknum(x)};
Cal.Date.CALC_DAYS=datenum(year_period(1),1,1):datenum(year_period(end),12,31);

for brw=1:length(Cal.brw)
    figure; set(gcf,'Tag',sprintf('Langley_R6_%s',Cal.brw_str{brw}));
    ha=gca;    idx=2; 
    
    % Configs: Operative
    try
       events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,Cal.brw_config_files{brw,1});    
       fprintf('\nBrewer %s: Operative Config.\n',Cal.brw_name{brw});
       displaytable(events_cfg_op.data(2:end,:),cellstr(datestr(events_cfg_op.data(1,:),1))',12,'.5g',events_cfg_op.legend(2:end));
    catch exception
       fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{brw});
    end
    mixed_brw_indv=cat(1,lang{brw}(:,[1 idx]),lang{brw}(:,[1 idx+1])); mixed_brw_indv=sortrows(mixed_brw_indv(~isnan(mixed_brw_indv(:,2)),:),1);
    dates=cell2mat(events_raw{brw}(:,2)); indx=dates>=lang{brw}(1,1) & dates<=lang{brw}(end,1); 
    [m_brw,s_brw,n_brw,std_brw]=grpstats(mixed_brw_indv,grp_def(mixed_brw_indv(:,1)),{'mean','sem','numel','std'});
    lmu=sortrows(m_brw,1); lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2)),1); lstd=sortrows(cat(2,m_brw(:,1),n_brw(:,2)),1);
  
    plot(lang{brw}(:,1),lang{brw}(:,idx),'b.',lang{brw}(:,1),lang{brw}(:,idx+1),'r.'); 
    boundedline(gca,lmu(:,1),lmu(:,2),lsem(:,2),'-k',lmu(:,1),lmu(:,2),lstd(:,2),'-k',...
                                              'alpha' ,'transparency', 0.3);
    ylabel('ETC (Brw method)');
    title(sprintf('Langley plot (%s - %s): %s. Config. Op.',...
              datestr(lang{brw}(1,1),28),datestr(lang{brw}(end,1),28),Cal.brw_name{brw})); 
  
    % Config from icf
    idx=group_time(Cal.Date.CALC_DAYS',events_cfg_op.data(1,:));
    stairs(gca,cat(2,events_cfg_op.data(1,logical(unique(idx))),lang{brw}(end,1)),cat(2,events_cfg_op.data(8,unique(idx)),events_cfg_op.data(8,end)),'-g','LineWidth',2);
      
    set(ha(1),'YLim',ylim_brw{brw}); linkprop(ha,'YLim');
    set(ha,'XLim',[Cal.Date.CALC_DAYS(1)-10 Cal.Date.CALC_DAYS(end)+10]); 
    datetick('x',12,'KeepLimits','keepTicks'); grid;    
end  