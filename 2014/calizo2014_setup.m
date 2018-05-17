%% SETUP
warning('off','MATLAB:xlsread:Mode');

%% General
% Habra que modificar a mano el directorio final: CODE, Are2011, SDK11 ...
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
Cal.file_save='calizo_2014.mat';
Cal.campaign='Izana 2014 (Canary Islands, Spain)';


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
Date.cal_year=2014;
Date.cal_month=04;

Cal.path_root=fullfile(Cal.path_root,'RBCC_E',num2str(Date.cal_year));

Date.CALC_DAYS=Date.day0:Date.dayend;
Date.range=[Date.CALC_DAYS(1),Date.CALC_DAYS(end)]; 
Date.BLIND_DAYS=Date.day0:Date.dayend;
Date.FINAL_DAYS=Date.day0:Date.dayend;

Date.N_Period=NaN;
Date.Period_start=NaN;
Date.Period.end=NaN;

Cal.Date=Date;

%% CALIBRATION INFO
Cal.Tsync=5;
Cal.brw=[157,183,185]; Cal.n_brw=length(Cal.brw);
Cal.brw_name={'IZO#157','IZO#183','IZO#185'};
Cal.brw_str=mmcellstr(sprintf('%03d|',Cal.brw));
Cal.brwM=[3,3,3];

Cal.brewer_ref=[1,2,3];
Cal.n_ref=[1,2,3];
Cal.no_maint=[0, 0, 0];
Cal.calibration_days={
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend  %157
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend  %183
   Date.day0:Date.dayend,Date.day0:Date.dayend,Date.day0:Date.dayend  %185
   %038:058,038:058,038:058                                            %221
   %090:140,[90:114 119 120 122:130 132:140],[90:114 119 120 122:130 132:140]  %145
   %315:323,315:323,315:323                                            %226
   %315:330,315:330,315:330                                            %227
                     };

% SL corr
Cal.sl_c      =[0,  0,  0];
Cal.sl_c_blind=[1,  1,  1];

%% Brewer configuration files. Eventos e Incidencias
icf_op=cell(1,length(Cal.n_brw)); icf_a=cell(1,length(Cal.n_brw));
events=cell(1,length(Cal.n_brw)); events_text=cell(1,length(Cal.n_brw)); events_raw=cell(1,length(Cal.n_brw));
incidences=cell(1,length(Cal.n_brw)); incidences_text=cell(1,length(Cal.n_brw)); incidences_raw=cell(1,length(Cal.n_brw));
for iz=1:Cal.n_brw
    icf_op{iz}=[]; icf_a{iz}=[];
 try   
    if iz<4 %reference
       icf_op{iz}=xlsread(fullfile(Cal.path_root,'..','configs',['icf',Cal.brw_str{iz},'.xls']),...
                         ['icf.',Cal.brw_str{iz}],'','basic');      
       icf_a{iz}=xlsread(fullfile(Cal.path_root,'..','configs',['icf',Cal.brw_str{iz},'.xls']),...
                         ['icf_a.',Cal.brw_str{iz}],'','basic');                     
       [events{iz},events_text{iz},events_raw{iz}]=xlsread(fullfile(Cal.path_root,'..','configs',['icf',Cal.brw_str{iz},'.xls']),...
                                                             ['Eventos.',Cal.brw_str{iz}],'','basic');
       [incidences{iz},incidences_text{iz},incidences_raw{iz}]=xlsread(fullfile(Cal.path_root,'..','configs',['icf',Cal.brw_str{iz},'.xls']),...
                                                             ['Incidencias.',Cal.brw_str{iz}],'','basic');
    else % not reference
      if exist(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']))
         icf_op{iz}=xlsread(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']),...
                           ['icf.',Cal.brw_str{iz}],'','basic');
         icf_a{iz}=xlsread(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']),...
                           ['icf_a.',Cal.brw_str{iz}],'','basic');                               
         [events{iz},events_text{iz},events_raw{iz}]=xlsread(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']),...
                                                               ['Eventos.',Cal.brw_str{iz}],'','basic');
          [incidences{iz},incidences_text{iz},incidences_raw{iz}]=xlsread(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},['icf',Cal.brw_str{iz},'.xls']),...
                                                 ['Incidencias.',Cal.brw_str{iz}],'','basic');         
      else
          continue
      end
    end   
    
    if size(icf_op{iz},1)==54  % ??
       cfg=icf_op{iz}(2:end-1,3:end); 
       save('config.cfg', 'cfg', '-ASCII','-double');
    else
       cfg=icf_op{iz}(1:end-1,3:end); 
       save('config.cfg', 'cfg', '-ASCII','-double');
    end
    tmp_file=sprintf('config%s.cfg',Cal.brw_str{iz});       
    copyfile('config.cfg',fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},tmp_file));
    delete('config.cfg');
  
    if size(icf_a{iz},1)==54 %??
       cfg=icf_a{iz}(2:end-1,3:end);
       save('config_a.cfg', 'cfg', '-ASCII','-double');
    else
       cfg=icf_a{iz}(1:end-1,3:end);
       save('config_a.cfg', 'cfg', '-ASCII','-double');
    end             
    tmp_file=sprintf('config%s_a.cfg',Cal.brw_str{iz});
    copyfile('config_a.cfg',fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},tmp_file));
    delete('config_a.cfg');
  
  catch exception
        fprintf('%s Brewer%s\n',exception.message,Cal.brw_name{iz});          
  end
  Cal.events{iz}=events{iz}(:,2:end);
  Cal.events_text{iz}=events_text{iz};
  Cal.events_raw{iz}=events_raw{iz};
  Cal.incidences_text{iz}=incidences_text{iz};
end

%% Operative is 1st, alternative is 2nd
  brw_config_files={
    '..\157\config157.cfg','..\157\config157_a.cfg','0350','0350';
    '..\183\config183.cfg','..\183\config183_a.cfg','0335','0335';
    '..\185\config185.cfg','..\185\config185_a.cfg','0305','0305';
    %'..\221\icf02714.221' ,'..\221\icf02714.221' ,'0325','0325';
    %'..\145\config145.cfg','..\145\config145_a.cfg' ,'0430','0430';
    %'..\226\icf29514.226' ,'..\226\icf31514.226' ,'0370','0290';
   % '..\227\icf29414.227' ,'..\227\icf31514.227' ,'0495','0512';
   };

  Cal.ETC_C={
             [0,0,0,0,0,0]          %157
             [0,0,0,0,-8,0]         %183
             [0,0,0,11,9,0]         %185
             %[0,0,0,0,0,0]          %221
             %[0,0,0,0,0,0]          %145
             %[0,0,0,15,15,0]        %226
             %[0,0,0,-6,-17,0]       %227
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
