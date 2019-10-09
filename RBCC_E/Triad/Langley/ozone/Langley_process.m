% options_pub=struct('outputDir',fullfile('..','..','html'),'showCode',false); 
% close all; publish(fullfile(pwd,'Langley157.m'),options_pub);

%%
% use matrix from load_data_2016_2019

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
%%


   nd0_=cat(1,brw_indv_{Cal.n_inst}(:,[1 2],cfgs),[brw_indv_{Cal.n_inst}(:,1,cfgs)+0.5,brw_indv_{Cal.n_inst}(:,[6],cfgs)]); nd0=sortrows(nd0_,1);
   nd3_=cat(1,brw_indv_{Cal.n_inst}(:,[1 3],cfgs),[brw_indv_{Cal.n_inst}(:,1,cfgs)+0.5,brw_indv_{Cal.n_inst}(:,7,cfgs)]); nd3=sortrows(nd3_,1);
   nd4_=cat(1,brw_indv_{Cal.n_inst}(:,[1 4],cfgs),[brw_indv_{Cal.n_inst}(:,1,cfgs)+0.5,brw_indv_{Cal.n_inst}(:,8,cfgs)]); nd4=sortrows(nd4_,1);
     %lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2:end)),1);
   brw_indv_filter{Cal.n_inst}=brw_indv_{Cal.n_inst};
  % save(Cal.file_save,'-APPEND','lgl_filt','brw_indv_filter');
   nd0_2=cat(1,brw_indv_{Cal.n_inst}(:,[1 2],cfgs+1),[brw_indv_{Cal.n_inst}(:,1,cfgs)+0.5,brw_indv_{Cal.n_inst}(:,[6],cfgs)]); nd0_2=sortrows(nd0_2,1);
   nd3_2=cat(1,brw_indv_{Cal.n_inst}(:,[1 3],cfgs+1),[brw_indv_{Cal.n_inst}(:,1,cfgs)+0.5,brw_indv_{Cal.n_inst}(:,[7],cfgs)]); nd3_2=sortrows(nd3_2,1);
   nd4_2=cat(1,brw_indv_{Cal.n_inst}(:,[1 4],cfgs+1),[brw_indv_{Cal.n_inst}(:,1,cfgs)+0.5,brw_indv_{Cal.n_inst}(:,[8],cfgs)]); nd4_2=sortrows(nd4_2,1);

  %  Salidas
           
    %t_days=array2table(days_lgl{Cal.n_inst}.data,'VariableNames',varname(days_lgl{Cal.n_inst}.labels),'Rownames',cellstr(datestr(days_lgl{Cal.n_inst}.data(:,1))));                  
    l_data_filter{Cal.n_inst}=array2table([nd0,nd3(:,2)-nd0(:,2),nd4(:,2)-nd0(:,2),nd0_2(:,2),nd3_2(:,2)-nd0_2(:,2),nd4_2(:,2)-nd0_2(:,2)],...
     'VariableNames',{'Diaj','ETC0','ND3','ND4','ETC0_cfg2','ND3_cfg2','ND4_cfg2'},...
     'Rownames',cellstr(datestr(nd0(:,1))));  
 
 
    writetable(l_data_filter{Cal.n_inst},'Langley_.xls','Sheet',strcat('filter_',Cal.brw_str{Cal.n_inst}))
    writetable(l_data_filter{Cal.n_inst},strcat('Langley_filter_',Cal.brw_str{Cal.n_inst}))

  
  

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
 

