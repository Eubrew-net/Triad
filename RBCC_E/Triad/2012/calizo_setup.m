%% SETUP

%% General
% Habrá que modificar a mano el directorio final: CODE, Are2011, SDK11 ...
if ispc
   Cal.path_root=fullfile(cell2mat(regexpi(pwd,'^[A-Z]:', 'match')),'CODE',''); 
else
   Cal.path_root=fullfile('~',cell2mat(regexpi(pwd,'^[A-Z]:', 'match')),'CODE','iberonesia'); 
end
path(genpath(fullfile(Cal.path_root,'matlab')),path);
Cal.file_save='calizo_2012.mat';
Cal.campaign='Izana (Spain) 2012';


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
Date.day0=1; Date.dayend=365;
Date.cal_year=2012; Date.cal_month=1;

Cal.path_root=fullfile(Cal.path_root,'iberonesia','RBCC_E',num2str(Date.cal_year)); 

Date.CALC_DAYS=Date.day0:Date.dayend;
Date.BLIND_DAYS=Date.day0:Date.dayend;
Date.FINAL_DAYS=Date.day0:Date.dayend;

Date.N_Period=NaN;
Date.Period_start=NaN;
Date.Period.end=NaN;

Cal.Date=Date;

%% CALIBRATION INFO
Cal.Tsync=5;
Cal.brw=[157,183,185,155,066,216]; Cal.n_brw=length(Cal.brw);
Cal.brw_name={'IZO#157','IZO#183','IZO#185','URU#155','AAB#066','K&Z#216'};
Cal.brw_str=mmcellstr(sprintf('%03d|',Cal.brw));
Cal.brwM=[3,3,3,2,4,3]; 

Cal.brewer_ref=[1,2,3];
Cal.n_ref=[1,2,3];
Cal.no_maint=1;
Cal.calibration_days={
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend         %157
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend         %183
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend         %185
   178:Date.dayend,178:Date.dayend,178:Date.dayend                           %155
   257:295,257:295,257:295                                                   %216
   340:347,340:341,342:347                                                   %066
                     };

% SL corr
Cal.sl_c=[0,  1,  0,  0, 0, 0];

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
  brw_config_files={% ios10 setting as stated in .\sdk2011\bfiles\185\Izana_cal2010.pdf
    '..\157\config157.cfg','..\157\config157_a.cfg','0350','0350';
    '..\183\config183.cfg','..\183\config183_a.cfg','0335','0335';
    '..\185\config185.cfg','..\185\config185_a.cfg','0210','0210';
    '..\155\icf15012.155' ,'..\155\config155.cfg'  ,'1540','1505';% icf27212.155 
    '..\066\icf21511.066' ,'..\066\icf20212.066'   ,'1840','1843'; 
    '..\216\icf31812.216' ,'..\216\icf34512.216'   ,'0295','0316';
    };

  Cal.ETC_C={
             [0,0,0,0,0,0]          %157
             [0,0,0,0,-8,0]         %183
             [0,0,0,11,9,0]         %185
             [0,0,0,0,0,0]          %155
             [0,0,0,18,-50,0]       %066
             [0,0,0,5,0,0]          %216
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