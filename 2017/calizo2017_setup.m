%% SETUP
warning('off','MATLAB:xlsread:Mode');
%thisYear=datestr(2017,'yyyy'); %J
thisYear='2017'; %J

%% General
% Habra que modificar a mano el directorio final: CODE, Are2011, SDK11 ...
if exist('alberto','file' )
       disp('alberto_setup')
       Cal.path_root=fullfile(cell2mat(regexpi(pwd,'^[A-Z]:', 'match')),'CODE','');
else
  if ispc
     Cal.path_root=fullfile(cell2mat(regexpi(pwd,'^[A-Z]:', 'match')),'CODE','iberonesia');
  else %J linux and macos
      [linuxStat,linuxUser]=system('whoami'); %J get the user name to determine its home
      linuxUser=strtrim(linuxUser); %J remove any extra spaces added by matlab
      if ismac
        Cal.path_root=fullfile(filesep,'Users',linuxUser,'CODE','rbcce.aemet.es','iberonesia')
      else    
      Cal.path_root=fullfile(filesep,'home',linuxUser,'CODE','iberonesia')
      end
  end
end

path(genpath(fullfile(Cal.path_root,'matlab')),path);
Cal.file_save=['calizo_' thisYear '.mat']; %J
Cal.campaign=['Izana ' thisYear ' (Canary Islands, Spain)']; %J


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
Date.cal_year=str2num(thisYear); %J

Date.day0=1;

if eomday(Date.cal_year,2)==29 %J
	Date.dayend=366;       %J
else                           %J
	Date.dayend=365;       %J
end                            %J



Date.cal_month=04; %J not used?
Date.year=2017;

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
Cal.brw=[157,183,185]; 
Cal.n_brw=length(Cal.brw);
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
    delete(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},tmp_file));
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
    delete(fullfile(Cal.path_root,'bfiles',Cal.brw_str{iz},tmp_file));
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
    %'..\157\icf05016.157','..\157\icf05116.157','0342','0342';
    '..\157\config157.cfg','..\157\config157_a.cfg','0342','0342';    
    '..\183\config183.cfg','..\183\config183_a.cfg','0395','0395';
    %'..\185\icf19016.185','..\185\icf11916.185','0312','0330';
    '..\185\config185.cfg','..\185\config185_a.cfg','0330','0330';
   };

  Cal.ETC_C={
             [0,0,0,0,0,0]          %157
             [0,0,0,0,-8,0]         %183
             [0,0,0,15,22,0]         %185
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

% for new hg graphics
maxf=@(x) x(end);
minf=@(x) x(1);
brewerdate=@(x,y) datenum(y,1,0)+x;
n1=@(x,y) unique(x(~isnan(x(:,y+1)),y+1));
%
Cal.maxf=maxf;
Cal.minf=minf;
Cal.brewerdate=brewerdate;
