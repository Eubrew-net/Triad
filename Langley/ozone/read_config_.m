
%'/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2018/bfiles/185/../185/config185.cfg'
file_setup='calizo_setup';
run(fullfile('..','..',file_setup));     % configuracion por defecto
op_cfg=cell(3,1);
chk_cfg=op_cfg;
events=cell(3,1);


Cal.Date.day0=datenum(2016,1,1);
Cal.Date.dayend=now; % to K&Z calibration
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend; Cal.Date=Date;

Cal.file_latex = fullfile('..','latex'); 
Cal.dir_tables = fullfile('..','latex','tables')
mkdir(Cal.file_latex);
Cal.dir_figs   = fullfile(Cal.file_latex,filesep(),'figures');
mkdir(Cal.dir_figs);
mkdir(Cal.dir_tables);


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
      op_cfg{i}=display_table(events_cfg_op.data(1:end,:),cellstr(datestr(events_cfg_op.data(1,:),1))',12,'.5g',events_cfg_op.legend(1:end));
      %Cal.n_inst=i;
      %events{i}=getevents(Cal,'grp','events');
    catch exception
      fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{i});
    end
     matrix2latex_config(events_cfg_op.data(2:end,:),fullfile(Cal.dir_tables,['Op_config_',Cal.brw_str{i},'.tex']),...
                     'rowlabels',str2var(events_cfg_op.legend(2:end)),'columnlabels',cellstr(datestr(events_cfg_op.data(1,:),1))',...
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
       Cal.n_inst=i;
       events{i}=getevents(Cal,'grp','events');
       alt_cfg{i}=display_table(events_cfg_chk.data(1:end,:),cellstr(datestr(events_cfg_chk.data(1,:),'mmmyy-dd'))',12,'.5g',events_cfg_chk.legend(1:end));
       
            matrix2latex_config(events_cfg_chk.data(2:end,:),fullfile(Cal.dir_tables,['ALT_config_',Cal.brw_str{i},'.tex']),...
                     'rowlabels',str2var(events_cfg_chk.legend(2:end)),'columnlabels',cellstr(datestr(events_cfg_chk.data(1,:),1))',...
                     'size','footnotesize');
       
     %%
     disp(alt_cfg{i})
     disp(events{i}.labels)
     
     try
          alt_cfg{i}.Properties.VariableNames=str2name(events{i}.labels)';
          
          
          disp(rows2vars(alt_cfg{i}(1,:)))
          op_cfg{i}.Properties.VariableNames=str2name(events{i}.labels)';
          disp(rows2vars(op_cfg{i}))
     catch
          warning('events and config not agree')
          events{i}.labels
          %events{i}.dates
          alt_cfg{i}(1,:)
     end
     
    catch exception
       fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{i});
    end
    
    
end
save   alt_cfg_16_19 alt_cfg
save   op_cfg_16_19  op_cfg
save   events_16_19 events
