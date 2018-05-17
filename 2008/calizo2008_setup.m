%% SETUP
warning('off','MATLAB:xlsread:Mode');

%% General
% Habrá que modificar a mano el directorio final: CODE, Are2011, SDK11 ...
if exist('alberto','file' )
       disp('alberto_setup')
       Cal.path_root=fullfile(cell2mat(regexpi(pwd,'^[A-Z]:', 'match')),'CODE','');
else
if ispc
   Cal.path_root=fullfile(cell2mat(regexpi(pwd,'^[A-Z]:', 'match')),'CODE','iberonesia');
else
   Cal.path_root=fullfile('~',cell2mat(regexpi(pwd,'^[A-Z]:', 'match')),'CODE','iberonesia');
end
end
path(genpath(fullfile(Cal.path_root,'matlab')),path);
Cal.file_save='calizo_2008.mat';
Cal.campaign='Izana (Spain)';


%% Station
Station.OSC=680;
Station.name='IZANA';
Station.lat=[];
Station.long=[];
Station.meanozo=[];

Cal.Station=Station;

%% Finalcalibration
Cal.FCal.ICF_FILE_INI='ICFXXXYY';
Cal.FCal.ICF_FILE_FIN='ICFXXXYY';
Cal.FCal.DCFFILE='DCFXXXYY';
Cal.FCal.LFFILE='LFXXXYY';

%%  configuration  date---> Default values
Date.day0=1;
Date.dayend=365;
Date.cal_year=2008;
Date.cal_month=09;

Cal.path_root=fullfile(Cal.path_root,'RBCC_E',num2str(Date.cal_year));

Date.CALC_DAYS=Date.day0:Date.dayend;
Date.BLIND_DAYS=Date.day0:Date.dayend;
Date.FINAL_DAYS=Date.day0:Date.dayend;

Date.N_Period=NaN;
Date.Period_start=NaN;
Date.Period.end=NaN;

Cal.Date=Date;

%% CALIBRATION INFO
Cal.Tsync=5;
Cal.brw=[157,183,185,008,033,065,145,193,196];
Cal.n_brw=length(Cal.brw);
Cal.brw_name={'IZO#157','IZO#183','IZO#185','EC_#008','SCO#033','EC_#145','K&Z#065','K&Z#193','K&Z#196'};
Cal.brw_str=mmcellstr(sprintf('%03d|',Cal.brw));
Cal.brwM=[3,3,3,2,2,4,3,2,3,3];

Cal.brewer_ref=[1,2,3];
Cal.n_ref=[1,2,3];
Cal.no_maint=1;
Cal.calibration_days={
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend  %157
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend  %183
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend  %185
   250:290,250:290, 250:290                                           %008
   283:Date.dayend, 284:286         , 287:290                         %033
   303:320        , 303:320         , 303:320                         %037
   303:321        , [303:317 319]   , [303:317 319]                   %214
   310:329        , [310:316 318:320 322:326] , [310:316 318:320 322:326]     %053
   316:329        , [316:317 320:329]         , [316:317 320:329]             %082
   320:329        , 320:329         , 320:329                         %202
                     };

% SL corr during calibration & during blind days
Cal.sl_c=[0,  0,  0, 0, 0, 1, 0, 0, 0, 0];
Cal.sl_c_blind=[1,  1,  1, 0, 0, 0, 0, 0, 0, 0];

%% Brewer configuration files. Eventos e Incidencias
icf_n={}; events_text={}; incidences_text={}; events_raw={};
for iz=1:Cal.n_brw
    %icf_n{iz}=[];
    if iz<4
      icf_n{iz}=xlsread(fullfile(Cal.path_root,'..','configs',['icf',Cal.brw_str{iz},'.xls']),...
                                    ['icf.',Cal.brw_str{iz}],'','basic');
      [events_n,events_text{iz},events_raw{iz}]=xlsread(fullfile(Cal.path_root,'..','configs',['icf',Cal.brw_str{iz},'.xls']),...
                                             ['Eventos.',Cal.brw_str{iz}],'','basic');
      [inc_n,incidences_text{iz},incidences_raw{iz}]=xlsread(fullfile(Cal.path_root,'..','configs',['icf',Cal.brw_str{iz},'.xls']),...
                                             ['Incidencias.',Cal.brw_str{iz}],'','basic');
    else
      if exist(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']))
        try
        icf_n{iz}=xlsread(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']),...
                                             ['icf.',Cal.brw_str{iz}],'','basic');
        [events_n,events_text{iz},events_raw{iz}]=xlsread(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']),...
                                             ['Eventos.',Cal.brw_str{iz}],'','basic');
        [inc_n,incidences_text{iz}]=xlsread(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']),...
                                             ['Incidencias.',Cal.brw_str{iz}],'','basic');
       catch exception
          fprintf('%s Brewer%s\n',exception.message,Cal.brw_name{iz});
        end
      else
          continue;
      end
    end

    if size(icf_n{iz},1)==54
       cfg=icf_n{iz}(2:end-1,3:end); save('config.cfg', 'cfg', '-ASCII','-double');
    else
       cfg=icf_n{iz}(1:end-1,3:end); save('config.cfg', 'cfg', '-ASCII','-double');
    end
    tmp_file=sprintf('config%s.cfg',Cal.brw_str{iz});
    copyfile('config.cfg',fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},tmp_file));
    delete('config.cfg');

    Cal.events{iz}=events_n(:,2:end);
    Cal.events_text{iz}=events_text{iz};
    Cal.events_raw{iz}=events_raw{iz};
    Cal.incidences_text{iz}=incidences_text{iz};
end

% Alternative configuration
for iz=1:3
    icfn=xlsread(fullfile(Cal.path_root,'..','configs',['icf',Cal.brw_str{iz},'.xls']),...
                   ['icf_a.',Cal.brw_str{iz}],'','basic');
    if size(icf_n{iz},1)==54
       cfg=icfn(2:end-1,3:end); save('config_a.cfg', 'cfg', '-ASCII','-double');
    else
       cfg=icfn(1:end-1,3:end); save('config_a.cfg', 'cfg', '-ASCII','-double');
    end
    tmp_file=sprintf('config%s_a.cfg',Cal.brw_str{iz});
    copyfile('config_a.cfg',fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},tmp_file));
    delete('config_a.cfg');
end

%%
  brw_config_files={
    '..\157\config157.cfg','..\157\config157.cfg','0350','0350';
    '..\183\config183.cfg','..\183\config183.cfg','0325','0325';
    '..\185\config185.cfg','..\185\config185.cfg','0305','0305';
    '..\008\icf20607.008','..\008\icf20607.008','1837','1837';
    '..\033\icf24207.033','..\033\icf24008.033','2130','2130';
    '..\145\icf33007.145','..\145\icf33007.145','0440','0440';
    '..\065\ICF24013.214' ,'..\214\ICF30313.214', '0357','0238';
    '..\193\ICF16111.053' ,'..\053\ICF32513.053', '1830','1845';
    '..\196\ICF18009.082' ,'..\082\ICF31613.082', '1555','1615';
   };

  Cal.ETC_C={
             [0,0,0,0,0,0]          %157
             [0,0,0,0,-8,0]         %183
             [0,0,0,11,9,0]         %185
             [0,0,0,0,0,0]          %008
             [0,0,0,0,0,0]          %033
             [0,0,0,0,0,0]          %145
             [0,0,0,0,0,0]          %065
             [0,0,0,0,0,0]          %193
             [0,0,0,0,0,0]          %196
            };

%%
pa=repmat(cellstr([Cal.path_root,filesep(),'bfiles']),Cal.n_brw,2);

pa=cellfun(@fullfile,pa,[mmcellstr(sprintf('%03d|',Cal.brw)),mmcellstr(sprintf('%03d|',Cal.brw))],'UniformOutput',0);
brw_config_files(:,1:2)=cellfun(@fullfile,pa,brw_config_files(:,1:2),'UniformOutput',0);
if isunix
 brw_config_files=strrep(brw_config_files,'\',filesep());
 %brw_config_files=cellfun(@upper,brw_config_files,'UniformOutput',0);
end
flag_config=cellfun(@exist,brw_config_files(:,1:2));
if ~all(all(flag_config))
    disp(flag_config);
    disp('Error config do not exist');
    brw_config_files(~flag_config)    
end

SL_OLD_REF=str2num(cat(1,brw_config_files{:,3}));
SL_NEW_REF=str2num(cat(1,brw_config_files{:,4}));
brw_config_files_old=brw_config_files(:,1);
brw_config_files_new=brw_config_files(:,2);

Cal.brw_config_files=brw_config_files;
Cal.brw_config_files_old=brw_config_files_old;
Cal.brw_config_files_new=brw_config_files_new;
Cal.SL_OLD_REF=str2num(cat(1,brw_config_files{:,3}));
Cal.SL_NEW_REF=str2num(cat(1,brw_config_files{:,4}));

%% Latex directories

%pa=repmat(cellstr([Cal.path_root,filesep(),'latex']),n_brw,1);
%pa=cellfun(@fullfile,pa,mmcellstr(sprintf('%03d|',Cal.brw)),'UniformOutput',0);
%Cal.dir_latex
%cellfun(@mkdir,pa)
%Cal.dir_figs=cellfun(@fullfile,pa,mmcellstr(sprintf('%03d_figures|',Cal.brw)),'UniformOutput',0);
%cellfun(@mkdir,Cal.Dir_figs)
