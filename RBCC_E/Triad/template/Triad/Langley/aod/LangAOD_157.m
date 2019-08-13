
Cal.n_inst=1;  cfgs=1;

airm_rang=[1.15 3.75];
grp_def=@(x) {year(x) weeknum(x)};

%% filters correction
ozone_lgl_dep{Cal.n_inst}=langley_filter_lvl1(ozone_lgl{Cal.n_inst},'plots',0,...
                             'airmass',airm_rang,'O3_hday',2);%,...
                            % 'AOD',fullfile('..','140101_141231_Izana.lev15'),...
                            % 'Cloud',fullfile('..','cloudScreening.txt'));

[brw_indv_aod_filters{Cal.n_inst} sptat_aod{Cal.n_inst}] = langley_analys_AOD_filter(ozone_lgl_dep,Cal.n_inst,Cal,...
                                                   'res_filt',1,'plot_flag',0);

%% No cloud data -> statistical filter 
for slit=1:5
    for filter=1:3
        % AM
       [ax,bx,cx,dx2]=outliers_bp(brw_indv_aod_filters{Cal.n_inst}(:,filter+1,slit,cfgs),1.5); 
        brw_indv_aod_filters{Cal.n_inst}(dx2,filter+1,slit,cfgs)=NaN;

        % PM
        [ax,bx,cx,dx2]=outliers_bp(brw_indv_aod_filters{Cal.n_inst}(:,filter+5,slit,cfgs),1.5); 
        brw_indv_aod_filters{Cal.n_inst}(dx2,filter+5,slit,cfgs)=NaN;
    end
end
 
% Write results to file
for cf=1:2
    aux=brw_indv_aod_filters{Cal.n_inst}(:,:,1,cf);
    for slit=2:5
        aux=cat(2,aux,brw_indv_aod_filters{Cal.n_inst}(:,2:end,slit,cf));
    end
    fid = fopen(sprintf('Brewer_LangAOD%d_%s_NDcorr_cfg%d.txt',Cal.Date.cal_year,Cal.brw_str{Cal.n_inst},cf), 'wt'); % Open for writing
    fprintf(fid, strcat('%%Date Slit#1_NDref(AM) Slit#1_ND#3(AM) Slit#1_ND#4(AM) Slit#1_slope(AM)',...
                               'Slit#1_NDref(PM) Slit#1_ND#3(PM) Slit#1_ND#4(PM) Slit#1_slope(PM)',...
                               'Slit#2_NDref(AM) Slit#2_ND#3(AM) Slit#2_ND#4(AM) Slit#2_slope(AM)',...
                               'Slit#2_NDref(PM) Slit#2_ND#3(PM) Slit#2_ND#4(PM) Slit#2_slope(PM)',...
                               'Slit#3_NDref(AM) Slit#3_ND#3(AM) Slit#3_ND#4(AM) Slit#3_slope(AM)',...
                               'Slit#3_NDref(PM) Slit#3_ND#3(PM) Slit#3_ND#4(PM) Slit#3_slope(PM)',...
                               'Slit#4_NDref(AM) Slit#4_ND#3(AM) Slit#4_ND#4(AM) Slit#4_slope(AM)',...
                               'Slit#4_NDref(PM) Slit#4_ND#3(PM) Slit#4_ND#4(PM) Slit#4_slope(PM)',...
                               'Slit#5_NDref(AM) Slit#5_ND#3(AM) Slit#5_ND#4(AM) Slit#5_slope(AM)',...
                               'Slit#5_NDref(PM) Slit#5_ND#3(PM) Slit#5_ND#4(PM) Slit#5_slope(PM)\n'));
    for i=1:size(aux,1)
        fprintf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n',aux(i,:));
    end
    fclose(fid);
end

%% Time Series
fprintf('\nNeutral Density Correction Factors (I0_corr): Brewer %s\r\n',Cal.brw_name{Cal.n_inst});

corr_slits=NaN*ones(2,5);  
for slit=1:5
    nd0_=cat(1,brw_indv_aod_filters{Cal.n_inst}(:,[1 2],slit,cfgs),brw_indv_aod_filters{Cal.n_inst}(:,[1 6],slit,cfgs)); nd0=sortrows(nd0_,1);
    nd3_=cat(1,brw_indv_aod_filters{Cal.n_inst}(:,[1 3],slit,cfgs),brw_indv_aod_filters{Cal.n_inst}(:,[1 7],slit,cfgs)); nd3=sortrows(nd3_,1);
    nd4_=cat(1,brw_indv_aod_filters{Cal.n_inst}(:,[1 4],slit,cfgs),brw_indv_aod_filters{Cal.n_inst}(:,[1 8],slit,cfgs)); nd4=sortrows(nd4_,1);
    
    figure; ha=tight_subplot(2,1,.048,.1,.075); hold all;
    axes(ha(1)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
    axes(ha(2)); set(gca,'box','on','YTickLabelMode','auto'); grid; hold on;
    
    % ND#3
    [m_brw,s_brw]=grpstats([nd0(:,1) nd3(:,2)-nd0(:,2)],grp_def(nd0(:,1)),{'mean','sem'});
    lmu=sortrows(m_brw,1); lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2:end)),1);   
    idx=lsem(:,2)==0 | isnan(lsem(:,2)); lsem(idx,:)=[];  lmu(idx,:)=[];  
    axes(ha(1)); [h p]=boundedline(gca,lmu(:,1),lmu(:,2),lsem(:,2),'-k','alpha' ,'transparency', 0.3);
    suptitle(sprintf('ND filters Corr. (Langley, %s - %s): %s, airmass range = [%.2f, %.2f]',...
             datestr(lmu(1,1),22),datestr(lmu(end,1),22),Cal.brw_name{Cal.n_inst},airm_rang)); 
    set(p,'Visible','off');
    h=hline(round(median(lmu(:,2))),'-r'); set(h,'LineWidth',2);
    legend(h,sprintf('Slit #%d: ND#3 Corr. = %d',slit,round(median(lmu(:,2)))));
    set(p,'Visible','on');
    corr_slits(1,slit)=median(lmu(:,2));        

    % ND#4: No measurements
    [m_brw,s_brw]=grpstats([nd0(:,1) nd4(:,2)-nd0(:,2)],grp_def(nd0(:,1)),{'mean','sem'});
    lmu=sortrows(m_brw,1); lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2:end)),1);   
    idx=lsem(:,2)==0 | isnan(lsem(:,2)); lsem(idx,:)=[];  lmu(idx,:)=[];  
    axes(ha(2)); [h p]=boundedline(gca,lmu(:,1),lmu(:,2),lsem(:,2),'-k','alpha' ,'transparency', 0.3);
    set(p,'Visible','off');
    h=hline(round(median(lmu(:,2))),'-r'); set(h,'LineWidth',2);
    legend(h,sprintf('Slit #%d: ND#4 Corr. = %d',slit,round(median(lmu(:,2)))));
    set(p,'Visible','on'); linkprop(ha,'XLim');  datetick('x',12,'KeepTicks');
    corr_slits(2,slit)=median(lmu(:,2));
end

%% O3W=[  0.00    0   0.00   -1.00    0.50    2.20   -1.70];
displaytable(corr_slits,{'Slit#1','Slit#2','Slit#3','Slit#4','Slit#5'},...
             12,'.5g',{'Filter corr.#3','Filter corr.#4'}); 
fprintf('\nMS9 corr. = %d\r\n',round(corr_slits*[0 -1.00    0.50    2.20   -1.70]'));

%% Procesamos los langleys individuales
close all
spectral_config=[0,0,0,0,0,0
                 0,0,0,0,0,0
                 0,0,0,0,0,0
                 0,0,0,0,0,0
                 0,0,0,0,0,0
                ];
ozone_lgl_dep{Cal.n_inst}=langley_filter_lvl1(ozone_lgl{Cal.n_inst},'plots',0,...
                             'airmass',airm_rang,'O3_hday',2,'F_corr_AOD',spectral_config);%,...
%                             'AOD',fullfile('..','140101_141231_Izana.lev15'));
%                              'date_range',datenum(2014,1,[200 290]));%,...
%                              'Cloud',fullfile('..','cloudScreening2013.txt'),...                              
%                              'lgl_days',1);

[brw_indv_aod{Cal.n_inst} stat_aod{Cal.n_inst}] = langley_analys_AOD(ozone_lgl_dep,Cal.n_inst,Cal,...
                                                   'res_filt',1,'plot_flag',0);
% Resumen de estadísticas
stats_summ=LangStats_summ(Cal,stat_aod{Cal.n_inst});

%% No cloud data -> statistical filter
for slit=1:5
    % AM
    [ax,bx,cx,dx2]=outliers_bp(brw_indv_aod{Cal.n_inst}(:,slit+1,cfgs),2.5);
    brw_indv_aod{Cal.n_inst}(dx2,slit+1,cfgs)=NaN;

    % PM
    [ax,bx,cx,dx2]=outliers_bp(brw_indv_aod{Cal.n_inst}(:,slit+6,cfgs),2.5);
    brw_indv_aod{Cal.n_inst}(dx2,slit+6,cfgs)=NaN;
end
% save(fullfile('..',Cal.file_save),'-APPEND','brw_indv_aod');

% Write results to file
for cf=1:2
    aux=brw_indv_aod{Cal.n_inst}(:,:,cf);
    fid = fopen(sprintf('Brewer_LangAOD%d_%s_cfg%d.txt',Cal.Date.cal_year,Cal.brw_str{Cal.n_inst},cf), 'wt'); % Open for writing
    fprintf(fid, strcat('%%Date Slit#1(AM) Slit#2(AM) Slit#3(AM) Slit#4(AM) Slit#5(AM)',...
                               'Slit#1(PM) Slit#2(PM) Slit#3(PM) Slit#4(PM) Slit#5(PM)',...
                               'Slope#1(AM) Slope#2(AM) Slope#3(AM) Slope#4(AM) Slope#5(AM)',...
                               'Slope#1(PM) Slope#2(PM) Slope#3(PM) Slope#4(PM) Slope#5(PM)\n'));
    for i=1:size(aux,1)
        fprintf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',aux(i,:));
    end
    fclose(fid);
end

%% Bounded plot
fprintf('\nIndividual I0: Brewer %s\r\n',Cal.brw_name{Cal.n_inst});

mixed_brw_indv=cell(1,5);

figure; set(gcf,'Tag',sprintf('Lag%s_bound',Cal.brw_str{Cal.n_inst}));
ha=tight_subplot(5,1,.048,.1,.075); hold all;
axes(ha(1)); set(ha(4),'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
axes(ha(2)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
axes(ha(3)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
axes(ha(4)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
axes(ha(5)); set(gca,'box','on','YTickLabelMode','auto'); grid; hold on;

for slit=1:5
    axes(ha(slit));
    plot(brw_indv_aod{Cal.n_inst}(:,1,cfgs),brw_indv_aod{Cal.n_inst}(:,slit+1,cfgs),'b.',...
         brw_indv_aod{Cal.n_inst}(:,1,cfgs),brw_indv_aod{Cal.n_inst}(:,slit+6,cfgs),'r.');

    mixed_brw_indv_=cat(1,brw_indv_aod{Cal.n_inst}(:,[1 slit+1],cfgs),brw_indv_aod{Cal.n_inst}(:,[1 slit+6],cfgs));
    mixed_brw_indv{slit}=sortrows(mixed_brw_indv_(~isnan(mixed_brw_indv_(:,2)),:),1);
    dates=cell2mat(events_raw{Cal.n_inst}(:,2)); indx=dates>=brw_indv_aod{Cal.n_inst}(1,1) & dates<=brw_indv_aod{Cal.n_inst}(end,1);

    [m_brw,s_brw,n_brw,std_brw]=grpstats(mixed_brw_indv{slit},grp_def(mixed_brw_indv{slit}(:,1)),{'mean','sem','numel','std'});
    lmu=sortrows(m_brw,1); lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2)),1); lstd=sortrows(cat(2,m_brw(:,1),n_brw(:,2)),1);
    boundedline(gca,lmu(:,1),lmu(:,2),lsem(:,2),'-k',lmu(:,1),lmu(:,2),lstd(:,2),'-k','alpha' ,'transparency', 0.3);
    set(gca,'YTicklabel',''); set(gca,'YTicklabelMode','manual','YTicklabel',get(ha(1),'YTick'));
    try
      vl=vline(dates(indx),'r-'); set(vl,'LineWidth',1);
    end
end
datetick('x',12,'keeplimits','keepticks');
title(ha(1),sprintf('Langley plot (%s - %s): %s, airmass range = [%.2f, %.2f]',...
      datestr(brw_indv_aod{Cal.n_inst}(1,1,cfgs),28),datestr(brw_indv_aod{Cal.n_inst}(end,1,cfgs),28),Cal.brw_name{Cal.n_inst},airm_rang));
ylabel(ha(3),'I0');
set(ha(1),'XLim',[brw_indv_aod{Cal.n_inst}(1,1,cfgs)-5 brw_indv_aod{Cal.n_inst}(end,1,cfgs)+5]); linkprop(ha,'XLim');

%% Brewer#157: ETC, Tabla por eventos + meses
data_=cell(1,5);
event_info=getevents(Cal,'grp','month+events');

fprintf('\r\nBrewer %s: ETC''s by slit (month + events average)\r\n',Cal.brw_name{Cal.n_inst});
for slit=1:5
    data_{slit}=meanperiods(mixed_brw_indv{slit}, event_info);
end
aux=NaN*ones(size(data_{1}.m,1),12);
aux(:,[2 4 6 8 10])=cell2mat(cellfun(@(x) x.m(:,2), data_,'UniformOutput',0));
aux(:,[3 5 7 9 11])=cell2mat(cellfun(@(x) x.std(:,2), data_,'UniformOutput',0));
aux(:,1)= data_{1}.m(:,1); aux(:,end)= data_{1}.N(:,end);

% Command-Print
displaytable(aux(:,2:end),{'Slit#1','#1std','Slit#2','#2std','Slit#3','#3std','Slit#4','#4std','Slit#5','#5std','N'},...
             8,'.1f',data_{1}.evnts);

%% Hist
% for slit=1:5
%     nd0_=cat(1,brw_indv_aod_filters{Cal.n_inst}(:,[1 2],slit,cfgs),brw_indv_aod_filters{Cal.n_inst}(:,[1 6],slit,cfgs)); nd0=sortrows(nd0_,1);
%     nd3_=cat(1,brw_indv_aod_filters{Cal.n_inst}(:,[1 3],slit,cfgs),brw_indv_aod_filters{Cal.n_inst}(:,[1 7],slit,cfgs)); nd3=sortrows(nd3_,1);
%     nd4_=cat(1,brw_indv_aod_filters{Cal.n_inst}(:,[1 4],slit,cfgs),brw_indv_aod_filters{Cal.n_inst}(:,[1 8],slit,cfgs)); nd4=sortrows(nd4_,1);
%
%     figure; %set(gcf,'Tag',sprintf('Lag%s_residuals',Cal.brw_str{Cal.n_inst}));
%     aux={nd0(:,2),nd3(:,2)};
%     nhist(aux,'box','smooth','samebins','ylabel','Pdf (Langley fit residuals)'); grid; box on;
%     set(findobj(gca,'Type','Line'),'Marker','None','LineWidth',2)
%     legendflex({'ND#0,1,2','ND#3'},'anchor', [2 6], 'buffer',[0 -20], ...
%                 'nrow',1,'fontsize',10,'box','on','xscale',.7,...
%                 'title',sprintf('Brewer %s: ND#3 Corr.=%d',Cal.brw_name{Cal.n_inst},...
%                         round(nanmedian(nd3(:,2))-nanmedian(nd0(:,2)))));
%     set(findobj(gcf,'Type','text'),'BackgroundColor','w');
% end

%% ploteo de los residuos para cada slit usando los dos DT
% Cal.n_inst=Cal.n_inst;
% for slit=5
%     % Op. Constants
%     brw_op{slit}=cell2mat(cellfun(@(x) x(~isnan(x(:,slit+3,1)),slit+3,1),stat_aod{Cal.n_inst}.r,'UniformOutput',0));
%     % Chk. Constants
%     brw_chk{slit}=cell2mat(cellfun(@(x) x(~isnan(x(:,slit+3,2)),slit+3,2),stat_aod{Cal.n_inst}.r,'UniformOutput',0));
% end
%
% for slit=5
%     figure; set(gcf,'Tag',sprintf('Lag%s_residuals',Cal.brw_str{Cal.n_inst}));
%     A_{slit}={brw_op{slit},brw_chk{slit}};
%     tex=nhist(A_{slit},'samebins','box','smooth','ylabel','Pdf (Langley fit residuals)');
% %     tex=nhist(A{slit}{1},'samebins','box');
% legend({'Brw Op. (DT=33)','Brw Chk. (DT=29)'},'location','NorthEast','FontSize',7);
% title(Cal.brw_name{Cal.n_inst}); set(findobj(gca,'Type','Line'),'Marker','None','LineWidth',1)
% end
