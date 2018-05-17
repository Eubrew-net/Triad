% options_pub.outputDir=fullfile(pwd,'html');   options_pub.showCode=false;
% publish(fullfile(pwd,'load_data_sergio2.m'),options_pub);

%%  Brewer setup
clear all;
file_setup='calizo2010_setup';
run(fullfile('..','..',file_setup));     % configuracion por defecto

Cal.Date.day0=datenum(2010,1,1);  Cal.Date.dayend=datenum(2010,12,31);
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend; Cal.Date=Date;

Cal.file_save  = 'triad_comp2010';
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
                               'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii},...
                               'hgflag',1,'fplot',1);
    % Imprimimos valores por eventos
    fprintf('SL means, Op. config (by events). Brewer %s\r\n',Cal.brw_name{ii}); Cal.n_inst=ii;
    event_info=getevents(Cal,'grp','events'); data_tab=meanperiods(sl_median_o{ii}(:,[1 2 4]), event_info);
    displaytable(cat(2,data_tab.m(:,2),data_tab.std(:,2),data_tab.m(:,3),data_tab.std(:,3),data_tab.N(:,1)),...
                 {'R6','std','R5','std','N'},15,'.4f',data_tab.evnts);
     % SL report plots
     ix=sort(findobj('Tag','SL_R6_report'));
     Width=15; Height=7;
     str=Cal.brw_str(Cal.n_ref);
     figure(ix);
     datetick('x',19,'keeplimits','keepticks');
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

     printfiles_report(gcf,Cal.dir_figs,'aux_pattern' ,str(ii),'Width',Width,'Height',Height,'LineMode','Scaled');
     snapnow; 
     close all
     
    [sl_mov_n{ii},sl_median_n{ii},sl_out_n{ii},R6_n{ii}]=sl_report_jday(ii,sl_cr,Cal.brw_str,...
                               'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii},...
                               'hgflag',1,'fplot',1);
    % Imprimimos valores por eventos
    fprintf('SL means, Chk. config (by events). Brewer %s\r\n',Cal.brw_name{ii}); Cal.n_inst=ii;
    event_info=getevents(Cal,'grp','events'); data_tab=meanperiods(sl_median_n{ii}(:,[1 2 4]), event_info);
    displaytable(cat(2,data_tab.m(:,2),data_tab.std(:,2),data_tab.m(:,3),data_tab.std(:,3),data_tab.N(:,1)),...
                 {'R6','std','R5','std','N'},15,'.4f',data_tab.evnts);
    % SL report plots
    ix=sort(findobj('Tag','SL_R6_report'));
    Width=15; Height=7;
    str=Cal.brw_str(Cal.n_ref);
    figure(ix);
    datetick('x',19,'keeplimits','keepticks');
    th=rotateticklabel(gca,45); set(th,'FontSize',9);
    set(findobj(gcf,'Type','text'),'FontSize',8); 
    set(gca,'XTicklabel','','FontSize',10); xlabel('');
    set(findobj(gca,'Marker','.'),'MarkerSize',5);
    set(findobj(gca,'Marker','o','-or','Marker','s'),'MarkerSize',4)
    set(findobj(gcf,'Tag','legend'),'FontSize',8,'Location','NorthWest','LineWidth',0.3);
    % Config from Operative icf
    ALT_config=Cal.brw_config_files{ii,2};
    events_cfg_al=getcfgs(Cal.Date.CALC_DAYS,ALT_config);
    idx=group_time(R6_n{ii}(:,1),events_cfg_al.data(1,:));
    stairs(gca,R6_n{ii}(logical(idx),1),events_cfg_al.data(17,idx),'-b','LineWidth',3);

    printfiles_report(gcf,Cal.dir_figs,'aux_pattern',str(ii),'Width',Width,'Height',Height,'LineMode','Scaled');
    snapnow;
    close all;         
disp(ii)
end
save(Cal.file_save,'-APPEND','sl_mov_o','sl_median_o','sl_out_o','R6_o',...
                             'sl_mov_n','sl_median_n','sl_out_n','R6_n');

% %% SL report plots
% ix=sort(findobj('Tag','SL_R6_report'));
% Width=15; Height=7;
% str=Cal.brw_str(Cal.n_ref);
% for i=1:length(ix)
%     figure(ix(i));
%     datetick('x',19,'keeplimits','keepticks');
%     th=rotateticklabel(gca,45); set(th,'FontSize',9);
%     set(findobj(gcf,'Type','text'),'FontSize',8); 
%     set(gca,'XTicklabel','','FontSize',10); xlabel('');
%     set(findobj(gca,'Marker','.'),'MarkerSize',5);
%     set(findobj(gca,'Marker','o','-or','Marker','s'),'MarkerSize',4)
%     if i~=2
%        set(findobj(gcf,'Tag','legend'),'FontSize',8,'Location','SouthWest','LineWidth',0.3);
%     else
%        set(findobj(gcf,'Tag','legend'),'FontSize',8,'Location','NorthWest','LineWidth',0.3);
%     end
%     % Config from Operative icf
%     OP_config=Cal.brw_config_files{i,1};
%     events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,OP_config);
%     idx=group_time(R6_n{i}(:,1),events_cfg_op.data(1,:));
%     stairs(gca,R6_n{i}(logical(idx),1),events_cfg_op.data(17,idx),'-b','LineWidth',3);
% 
%     printfiles_report(gcf,Cal.dir_figs,'aux_pattern',str(i),'Width',Width,'Height',Height,'LineMode','Scaled');
% end

% close all;
% for ii=Cal.n_ref
%     sl_mov_o{ii}={}; sl_median_o{ii}={}; sl_out_o{ii}={}; R6_o{ii}={};
%     sl_mov_n{ii}={}; sl_median_n{ii}={}; sl_out_n{ii}={}; R6_n{ii}={};
% % old instrumental constants
%     [sl_mov_o{ii},sl_median_o{ii},sl_out_o{ii},R6_o{ii}]=sl_report_jday(ii,sl,Cal.brw_str,...
%                                'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii},...
%                                'hgflag',1,'fplot',1);
%  % Imprimimos valores por eventos
%     fprintf('SL means, Op. config (by events). Brewer %s\r\n',Cal.brw_name{ii}); Cal.n_inst=ii;
%     event_info=getevents(Cal,'grp','events'); data_tab=meanperiods(sl_median_o{ii}(:,[1 2 4]), event_info);
%     displaytable(cat(2,data_tab.m(:,2),data_tab.std(:,2),data_tab.m(:,3),data_tab.std(:,3),data_tab.N(:,1)),...
%                  {'R6','std','R5','std','N'},15,'.4f',data_tab.evnts);
% 
% % new instrumental constants
%     [sl_mov_n{ii},sl_median_n{ii},sl_out_n{ii},R6_n{ii}]=sl_report_jday(ii,sl_cr,Cal.brw_str,...
%                                'outlier_flag',1,'diaj_flag',0,'events_raw',events_raw{ii},...
%                                'hgflag',1,'fplot',0);
% % Imprimimos valores por eventos
%     fprintf('SL means, Chk. config (by events). Brewer %s\r\n',Cal.brw_name{ii}); Cal.n_inst=ii;
%     event_info=getevents(Cal,'grp','events'); data_tab=meanperiods(sl_median_n{ii}(:,[1 2 4]), event_info);
%     displaytable(cat(2,data_tab.m(:,2),data_tab.std(:,2),data_tab.m(:,3),data_tab.std(:,3),data_tab.N(:,1)),...
%                  {'R6','std','R5','std','N'},15,'.4f',data_tab.evnts);
% 
% end
% save(Cal.file_save,'-APPEND','sl_mov_o','sl_median_o','sl_out_o','R6_o',...
%                              'sl_mov_n','sl_median_n','sl_out_n','R6_n');
% 
% %% SL report plots
% ix=sort(findobj('Tag','SL_R6_report')); Width=15; Height=7; str=Cal.brw_str(Cal.n_ref);
% for i=1:length(ix)
%     figure(ix(i));
%     datetick('x',19,'keeplimits','keepticks'); th=rotateticklabel(gca,45); set(th,'FontSize',9);
%     set(findobj(gcf,'Type','text'),'FontSize',8); set(gca,'XTicklabel','','FontSize',10); xlabel('');
%     set(findobj(gca,'Marker','.'),'MarkerSize',5); set(findobj(gca,'Marker','o','-or','Marker','s'),'MarkerSize',4)
%     if i~=2
%        set(findobj(gcf,'Tag','legend'),'FontSize',8,'Location','SouthWest','LineWidth',0.3);
%     else
%        set(findobj(gcf,'Tag','legend'),'FontSize',8,'Location','NorthWest','LineWidth',0.3);
%     end
%     % Config from Operative icf
%     OP_config=Cal.brw_config_files{i,1};
%     events_cfg_op=getcfgs(Cal.Date.CALC_DAYS,OP_config);
%     idx=group_time(R6_n{i}(:,1),events_cfg_op.data(1,:));
%     stairs(gca,R6_n{i}(logical(idx),1),events_cfg_op.data(17,idx),'-b','LineWidth',3);
% 
%     printfiles_report(gcf,Cal.dir_figs,'aux_pattern',str(i),'Width',Width,'Height',Height,'LineMode','Scaled');
% end

%% Creating Summaries
close all
[A,ETC,SL_B,SL_R,F_corr,SL_corr_flag,cfg]=read_cal_config_new(config,Cal,{sl_median_o,sl_median_n});

% Data recalculation for summaries  and individual observations
for i=Cal.n_ref
    cal{i}={}; summary_orig{i}={}; summary_orig_old{i}={};
    [cal{i},summary_orig{i},summary_orig_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,...
                                                      'flag_sl',0,'plot_sl',1,'flag_sl_corr',SL_corr_flag);
    % filter correction
    [summary_old{i} summary{i}]=filter_corr(summary_orig,summary_orig_old,i,A,F_corr{i});
end
write_summary(Cal.brw(1:3),Cal.Date.cal_year,summary_old,summary,SL_R,SL_B);
save(Cal.file_save,'-APPEND','A','ETC','F_corr','SL_B','SL_R','SL_corr_flag','cfg',...
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

% Introducido por Sergio para cargarme los outliers

% detectamos cuantos NaN tenemos:
tt=size(ref);
ref_new=zeros(tt(1,1),tt(1,2));
ref_new(:,1)=ref(:,1);

for d=1:1:tt(1,1)
    % detectamos cuantos NaN tenemos en cada fila de datos.
    dd=isnan(ref(d,:));
    ddd=dd(1,2)+dd(1,3)+dd(1,4);
    if ddd==0 % No hay NaN
        diff12=abs(ref(d,2)-ref(d,3));
        if diff12<=3.999999
            diff23=abs(ref(d,3)-ref(d,4));
            if diff23<=4
                ref_new(d,2)=ref(d,2);
                ref_new(d,3)=ref(d,3);
                ref_new(d,4)=ref(d,4);
            else
                ref_new(d,2)=ref(d,2);
                ref_new(d,3)=ref(d,3);
                ref_new(d,4)=NaN;
            end
        end
        
        if diff12>=4
            diff13=abs(ref(d,2)-ref(d,4));
            diff23=abs(ref(d,3)-ref(d,4));
            %Caso1: Brewer 1 y 3 calibrados y brewer 2 descalibrado.
            if diff13<=4 & diff23>=4
                ref_new(d,2)=ref(d,2);
                ref_new(d,3)=NaN;
                ref_new(d,4)=ref(d,4);
            end
            %Caso2: Brewer 2 y 3 calibrados y brewer 1 descalibrado.
            if diff23<=4 & diff13>=4
                ref_new(d,2)=NaN;
                ref_new(d,3)=ref(d,3);
                ref_new(d,4)=ref(d,4);
            end
        end
    end
   
    if ddd==1
        indice=0;
        for jj=2:1:4
            if dd(1,jj)==1
                indice=jj;
            end
        end
        if indice==2
            diff=abs(ref(d,3)-ref(d,4));
            if diff<=4
                ref_new(d,2)=NaN;
                ref_new(d,3)=ref(d,3);
                ref_new(d,4)=ref(d,4);
            else
                ref_new(d,2)=NaN;
                ref_new(d,3)=NaN;
                ref_new(d,4)=NaN;
            end
        end
        if indice==3
            diff=abs(ref(d,2)-ref(d,4));
            if diff<=4
                ref_new(d,2)=ref(d,2);
                ref_new(d,3)=NaN;
                ref_new(d,4)=ref(d,4);
            else
                ref_new(d,2)=NaN;
                ref_new(d,3)=NaN;
                ref_new(d,4)=NaN;
            end
        end
        if indice==4
            diff=abs(ref(d,2)-ref(d,3));
            if diff<=4
                ref_new(d,2)=ref(d,2);
                ref_new(d,3)=ref(d,3);
                ref_new(d,4)=NaN;
            else
                ref_new(d,2)=NaN;
                ref_new(d,3)=NaN;
                ref_new(d,4)=NaN;
            end
        end
    end
    if ddd==2
        ref_new(d,2)=NaN;
        ref_new(d,3)=NaN;
        ref_new(d,4)=NaN;
    end
    if ddd==3
        ref_new(d,2)=NaN;
        ref_new(d,3)=NaN;
        ref_new(d,4)=NaN;
    end            
end


figure; set(gcf,'Tag','Ozone No outliers');
ploty(ref_new,'.'); grid
ylim([200 400]);
legend(gca,Cal.brw_name{Cal.n_ref},'Location','NorthEast','Orientation','Horizontal');
datetick('x',6,'Keeplimits','Keepticks'); title('DS Ozone'); ylabel('Ozone (DU)');

clear d;clear dd;clear ddd; clear diff12; clear diff23; clear diff13;

% fid = fopen('rosa_data.txt', 'wt'); % Open for writing
% fprintf(fid, '%%Date Ozone(157) Ozone(183) Ozone(185)\n');
% for i=1:size(ref,1)
%     fprintf(fid, '%f %5.1f %5.1f %5.1f\n', ref(i,:));
% end
% fclose(fid);

%% Comparacion Operativa. (Todo el año)
reference_brw=[1 2 3]; analyzed_brewer=[1 2 3];
osc_interval=[400,700,1000,1200];
Cal.analyzed_brewer=analyzed_brewer; Cal.sl_c=[1,1,1];

Cal.brw=[157 183 185];
[ref,ratio_ref]=join_summary(Cal,summary_old,reference_brw,analyzed_brewer,3);
% 
% % osc table
% [ratio_osc_table,osc_matrix,osc_stats]=osc_table(Cal,ratio_ref,osc_interval);
% set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% 
% % ratio_ref_plots


% depuramos los erróres que hay en la variable "ratio_ref".
d=size(ratio_ref);
contador=0;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        contador=contador+1;
    end
end
Validos=NaN(contador,5);
rechazados=NaN((d(1,1)-contador),5);
cont_val=1;
cont_rech=1;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if  diff<1.2
        Validos(cont_val,:)=ratio_ref(dd,:);
        cont_val=cont_val+1;
    else
        rechazados(cont_rech,:)=ratio_ref(dd,:);
        cont_rech=cont_rech+1;
    end
end

[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,Validos,'plot_smooth',1);
fi=fopen('Todo_año_operativa_Validos.txt','w');
for kk=1:1:contador
    fprintf(fi,'%f %f %f %f %f/n',Validos(kk,1),Validos(kk,2),Validos(kk,3),Validos(kk,4),Validos(kk,5));
end
fclose(fi);
fi=fopen('Todo_año_operativa_rechazados.txt','w');
for kk=1:1:(d(1,1)-contador)
    fprintf(fi,'%f %f %f %f %f/n',rechazados(kk,1),rechazados(kk,2),rechazados(kk,3),rechazados(kk,4),rechazados(kk,5));
end
fclose(fi);

clear Validos
%% BEFORE AROSA
[ref,ratio_ref]=join_summary(Cal,summary_old,reference_brw,analyzed_brewer,5,'date_range',datenum(2010,1,[1,200]));
d=size(ratio_ref);
contador=0;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        contador=contador+1;
    end
end
Validos=NaN(contador,5);
rechazados=NaN((d(1,1)-contador),5);
cont_val=1;
cont_rech=1;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        Validos(cont_val,:)=ratio_ref(dd,:);
        cont_val=cont_val+1;
    else
        rechazados(cont_rech,:)=ratio_ref(dd,:);
        cont_rech=cont_rech+1;
    end
end
% ratio_ref_plots
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,Validos,'plot_smooth',1);
fi=fopen('Before_arosa_operativa_Validos.txt','w');
for kk=1:1:contador
    fprintf(fi,'%f %f %f %f %f/n',Validos(kk,1),Validos(kk,2),Validos(kk,3),Validos(kk,4),Validos(kk,5));
end
fclose(fi);
fi=fopen('Before_arosa_operativa_rechazados.txt','w');
for kk=1:1:(d(1,1)-contador)
    fprintf(fi,'%f %f %f %f %f/n',rechazados(kk,1),rechazados(kk,2),rechazados(kk,3),rechazados(kk,4),rechazados(kk,5));
end
fclose(fi);
clear Validos

%%  DURING ARENOSILLO

[ref,ratio_ref]=join_summary(Cal,summary_old,1:2,1:2,5,'date_range',datenum(2010,1,[201,213]));
Cal.analyzed_brewer=[1,2];
d=size(ratio_ref)
contador=0
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3));
    if diff<1.2
        contador=contador+1;
    end
end
Validos=NaN(contador,4);
rechazados=NaN((d(1,1)-contador),4);
cont_val=1;
cont_rech=1;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3));
    if diff<1.2
        Validos(cont_val,:)=ratio_ref(dd,:);
        cont_val=cont_val+1;
    else
        rechazados(cont_rech,:)=ratio_ref(dd,:);
        cont_rech=cont_rech+1;
    end
end
% ratio_ref_plots
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,Validos,'plot_smooth',1);
fi=fopen('During_arosa_operativa_Validos.txt','w');
for kk=1:1:contador
    fprintf(fi,'%f %f %f %f %f/n',Validos(kk,1),Validos(kk,2),Validos(kk,3),Validos(kk,4));
end
fclose(fi);
fi=fopen('During_arosa_operativa_rechazados.txt','w');
for kk=1:1:(d(1,1)-contador)
    fprintf(fi,'%f %f %f %f %f/n',rechazados(kk,1),rechazados(kk,2),rechazados(kk,3),rechazados(kk,4));
end
fclose(fi);
%% AFTER AROSA
[ref,ratio_ref]=join_summary(Cal,summary_old,1:3,1:3,5,'date_range',datenum(2010,1,[214,365]));
d=size(ratio_ref);
contador=0;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        contador=contador+1;
    end
end
Validos=NaN(contador,5);
rechazados=NaN((d(1,1)-contador),5);
cont_val=1;
cont_rech=1;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        Validos(cont_val,:)=ratio_ref(dd,:);
        cont_val=cont_val+1;
    else
        rechazados(cont_rech,:)=ratio_ref(dd,:);
        cont_rech=cont_rech+1;
    end
end
% ratio_ref_plots
Cal.analyzed_brewer=[1,2,3];
Cal.reference_brewer=[1,2,3];
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,Validos,'plot_smooth',1);
fi=fopen('After_arosa_operativa_Validos.txt','w');
for kk=1:1:contador
    fprintf(fi,'%f %f %f %f %f/n',Validos(kk,1),Validos(kk,2),Validos(kk,3),Validos(kk,4),Validos(kk,5));
end
fclose(fi);
fi=fopen('After_arosa_operativa_rechazados.txt','w');
for kk=1:1:(d(1,1)-contador)
    fprintf(fi,'%f %f %f %f %f/n',rechazados(kk,1),rechazados(kk,2),rechazados(kk,3),rechazados(kk,4),rechazados(kk,5));
end
fclose(fi);
%% Comparacion Alternativa
% close all

reference_brw=[1 2 3]; analyzed_brewer=[1 2 3];
osc_interval=[400,700,1000,1200];
Cal.analyzed_brewer=analyzed_brewer; Cal.sl_c=[0,1,0];

[ref,ratio_ref]=join_summary(Cal,summary,reference_brw,analyzed_brewer,10);
% 
% % osc table
% [ratio_osc_table,osc_matrix,osc_stats]=osc_table(Cal,ratio_ref,osc_interval);
% set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
% 
% % ratio_ref_plots
d=size(ratio_ref);
contador=0;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        contador=contador+1;
    end
end
Validos=NaN(contador,5);
rechazados=NaN((d(1,1)-contador),5);
cont_val=1;
cont_rech=1;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        Validos(cont_val,:)=ratio_ref(dd,:);
        cont_val=cont_val+1;
    else
        rechazados(cont_rech,:)=ratio_ref(dd,:);
        cont_rech=cont_rech+1;
    end
end
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,Validos,'plot_smooth',1);
fi=fopen('Todo_año_alternativa_Validos.txt','w');
for kk=1:1:contador
    fprintf(fi,'%f %f %f %f %f/n',Validos(kk,1),Validos(kk,2),Validos(kk,3),Validos(kk,4),Validos(kk,5));
end
fclose(fi);
fi=fopen('Todo_año_alternativa_rechazados.txt','w');
for kk=1:1:(d(1,1)-contador)
    fprintf(fi,'%f %f %f %f %f/n',rechazados(kk,1),rechazados(kk,2),rechazados(kk,3),rechazados(kk,4),rechazados(kk,5));
end
fclose(fi);
%% BEFORE ARENOSILLO
[ref,ratio_ref]=join_summary(Cal,summary_old,reference_brw,analyzed_brewer,5,'date_range',datenum(2010,1,[1,200]));
d=size(ratio_ref);
contador=0;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        contador=contador+1;
    end
end
Validos=NaN(contador,5);
rechazados=NaN((d(1,1)-contador),5);
cont_val=1;
cont_rech=1;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        Validos(cont_val,:)=ratio_ref(dd,:);
        cont_val=cont_val+1;
    else
        rechazados(cont_rech,:)=ratio_ref(dd,:);
        cont_rech=cont_rech+1;
    end
end
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,Validos,'plot_smooth',1);
fi=fopen('Before_arosa_alternativa_Validos.txt','w');
for kk=1:1:contador
    fprintf(fi,'%f %f %f %f %f/n',Validos(kk,1),Validos(kk,2),Validos(kk,3),Validos(kk,4),Validos(kk,5));
end
fclose(fi);
fi=fopen('Before_arosa_alternativa_rechazados.txt','w');
for kk=1:1:(d(1,1)-contador)
    fprintf(fi,'%f %f %f %f %f/n',rechazados(kk,1),rechazados(kk,2),rechazados(kk,3),rechazados(kk,4),rechazados(kk,5));
end
fclose(fi);

%%  DURING ARENOSILLO

[ref,ratio_ref]=join_summary(Cal,summary_old,1:2,1:2,5,'date_range',datenum(2010,1,[201,213]));
% ratio_ref_plots
Cal.analyzed_brewer=[1,2];
d=size(ratio_ref);
contador=0;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3));
    if diff<1.2
        contador=contador+1;
    end
end
Validos=NaN(contador,4);
rechazados=NaN((d(1,1)-contador),4);
cont_val=1;
cont_rech=1;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3));
    if diff<1.2
        Validos(cont_val,:)=ratio_ref(dd,:);
        cont_val=cont_val+1;
    else
        rechazados(cont_rech,:)=ratio_ref(dd,:);
        cont_rech=cont_rech+1;
    end
end
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,Validos,'plot_smooth',1);
fi=fopen('During_arosa_alternativa_Validos.txt','w');
for kk=1:1:contador
    fprintf(fi,'%f %f %f %f %f/n',Validos(kk,1),Validos(kk,2),Validos(kk,3),Validos(kk,4));
end
fclose(fi);
fi=fopen('During_arosa_alternativa_rechazados.txt','w');
for kk=1:1:(d(1,1)-contador)
    fprintf(fi,'%f %f %f %f %f/n',rechazados(kk,1),rechazados(kk,2),rechazados(kk,3),rechazados(kk,4));
end
fclose(fi);

%% AFTER ARENOSILLO
[ref,ratio_ref]=join_summary(Cal,summary_old,1:3,1:3,5,'date_range',datenum(2010,1,[214,365]));
% ratio_ref_plots
Cal.analyzed_brewer=[1,2,3];
Cal.reference_brewer=[1,2,3];
d=size(ratio_ref);
contador=0;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        contador=contador+1;
    end
end
Validos=NaN(contador,5);
rechazados=NaN((d(1,1)-contador),5);
cont_val=1;
cont_rech=1;
for dd=1:1:d(1,1)
    diff=abs(ratio_ref(dd,2))+abs(ratio_ref(dd,3))+abs(ratio_ref(dd,4));
    if diff<1.2
        Validos(cont_val,:)=ratio_ref(dd,:);
        cont_val=cont_val+1;
    else
        rechazados(cont_rech,:)=ratio_ref(dd,:);
        cont_rech=cont_rech+1;
    end
end
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
[f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,Validos,'plot_smooth',1);
fi=fopen('After_arenosillo_alternativa_validos.txt','w');
for kk=1:1:contador
    fprintf(fi,'%f %f %f %f %f/n',Validos(kk,1),Validos(kk,2),Validos(kk,3),Validos(kk,4),Validos(kk,5));
end
fclose(fi);
fi=fopen('After_arenosillo_alternativa_rechazados.txt','w');
for kk=1:1:(d(1,1)-contador)
    fprintf(fi,'%f %f %f %f %f/n',rechazados(kk,1),rechazados(kk,2),rechazados(kk,3),rechazados(kk,4),rechazados(kk,5));
end
fclose(fi);

%% Comparativa método triada Cánada

% Daily model for midday
% O3= Ai*Ii + B*(t-to) +C(t-to)^2
% to=solar noon


sum_day=[];
for dia=1:365
r_date=datenum(2010,0,0)+dia;
noon_m=solar_noon(r_date, -16.499);
r_date=r_date+noon_m/24/60;
diag=[];
for j=1:3

  if dia>200 & dia<214 & j==3
      
  cdata=1;
  else 
    jd=find(diaj(summary{j}(:,1))==dia);
    idx=zeros(size(jd,1),3);
    idx(:,j)=1;
    diag_=[summary{j}(jd,6),idx,summary{j}(jd,1)-r_date,(summary{j}(jd,1)-r_date).^2,(summary{j}(jd,1)-r_date).^3];
    diag=[diag;diag_];
  end
end
[B,BINT,R,RINT,STATS] = regress(diag(:,1),diag(:,2:end));
sum_day=[sum_day;[dia,B',STATS]];

end
sum_day(sum_day==0)=nan;
sum_day(abs(sum_day)>1000)=nan;

r1=mean(sum_day(:,2:4),2);
figure
plot(sum_day(:,1),100*matdiv(matadd(sum_day(:,2:4),-r1),r1))
ylim([-1.25 1.25])
save workspace