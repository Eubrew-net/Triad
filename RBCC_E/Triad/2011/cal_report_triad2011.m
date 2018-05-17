% options_pub.outputDir=fullfile('..','html','2009'); publish(fullfile(pwd,'cal_report_triad.m'),options_pub);

%% Brewer Evaluation
clear all; 
path(genpath(fullfile(pwd,'../matlab')),path);
file_setup='join_setup'; 

eval(file_setup);     % configuracion por defecto
Cal.n_inst=find(Cal.brw==185);
Cal.file_latex=fullfile('.','latex',Cal.brw_str{Cal.n_inst});
Cal.dir_figs=fullfile('latex',filesep(),Cal.brw_str{Cal.n_inst},...
                              filesep(),[Cal.brw_str{Cal.n_inst},'_figures'],filesep());
mkdir(Cal.dir_figs)
Cal.file_save='Triad2011';
try
 save(Cal.file_save,'-Append','Cal'); %sobreescribimos la configuracion guardada.   
 load(Cal.file_save); 
catch
    disp('clean');
    save(Cal.file_save);
end
disp(Cal.n_inst)

%% Brewer setup

Station.OSC=680;
Station.name='';
Station.lat=28;
Station.long=-18.5;
Station.meanozo=320;
% Usamos para la calibración días 15/Abr a 15/Jun (R6 estable y antes de cambio en calibración del 157)

Date.day0=-65; Date.dayend=diaj(now);
Date.CALC_DAYS=Date.day0:Date.dayend;
Date.BLIND_DAYS=204:206;
Date.FINAL_DAYS=Date.day0:Date.dayend;

latexcmd(fullfile(Cal.file_latex,['cal_setup_',Cal.brw_str{Cal.n_inst}]),...
                    '\CALINI',Date.CALC_DAYS(1),'\CALEND',Date.CALC_DAYS(end),...
                    '\slref',Cal.SL_OLD_REF(Cal.n_inst),'\slrefNEW',Cal.SL_NEW_REF(Cal.n_inst),...
                    '\BLINDINI',Date.BLIND_DAYS(1),'\BLINDEND',Date.BLIND_DAYS(end),...      
                    '\FINALINI',Date.FINAL_DAYS(1),'\FINALEND',Date.FINAL_DAYS(end),...
                    '\caldays',length(Date.FINAL_DAYS),'\Tsync',Cal.Tsync,...
                    '\brwname',Cal.brw_name(Cal.n_inst),'\brwref',Cal.brw_name(Cal.n_ref(2)),...
                    '\BRWSTATION',Station.name,'\STATIONOSC',Station.OSC,...
                    '\DCFFILE',Cal.FCal.DCFFILE,'\LFFILE',Cal.FCal.LFFILE);

Cal.Date=Date;
save(Cal.file_save,'-Append','Cal'); 

%% READ Brewer Summaries
% variable init
dsum=cell(Cal.n_brw,1);ozone_s={};ozone_ds={};
ozone_sum={};ds={};config={};sl_cr={}; dsum={};
ozone_sum={};config={};ozone_ds={}; com={};ds={}; 
missing=cell(Cal.n_brw);log=cell(Cal.n_brw);
missing=cell(Cal.n_brw);log=cell(Cal.n_brw);

 for i=1:Cal.n_brw
     dsum{i}={};    ozone_sum{i}={};  ozone_ds{i}={};
     config{i}={};    com{i}={};     ds{i}={};
     sl{i}={};     sl_cr{i}={};     hg{i}={};     bhg{i}={};
    missing{i}=[] ;log{i}={};
    
    [ozone,log_,missing_]=read_bdata(i,Cal);
    
    dsum{i}=ozone.dsum;
    ozone_sum{i}=ozone.ozone_sum;
    config{i}=ozone.config;
    ozone_ds{i}=ozone.ozone_ds;
    ozone_raw{i}=ozone.ozone_ds;
    hg{i}=ozone.hg;
    bhg{i}=ozone.bhg;
    sl{i}=ozone.sl; %first calibration/ bfiles
    sl_cr{i}=ozone.sl_cr; %recalculated with 2º configuration
    log{i}=cat(1,log_{:});
    missing{i}=missing_';

    disp(log{i});
 end

save(Cal.file_save,'-APPEND','ozone_raw','dsum','ozone_sum','ozone_ds','config','sl','sl_cr','hg','log','missing');
matrix2latex(log{Cal.n_inst},fullfile(Cal.file_latex,['tabla_fileprocess_',Cal.brw_str{Cal.n_inst},'.tex']));

%% HG Script en pruebas para desarrollar la funcion
% input celda hg,
% fecha matlab, hora , minituos ,segundos, correlacion
  % paso calculado, paso real,intensidad, temperatura,diff  pasos
  % cal-step , flag (?)
  %plot((hg1(:,1)),hg1(:,6)-mean(hg1(:,6)),'s')
  %hold on;
 
  hg_leg={'fecha','hora','minituos','segundos','correlacion',...     
  'step calc', 'step ','hg int','temp','step diff','hg flag'};
  % cal-step , flag (?)
 
  %nos interesa, paso, paso real, diff paso intensidad temperatura
  col_p=[6,7,8,10,5,9]
 
for i=1:size(hg,2)
  hg1=cell2mat(hg{i}(1:end));
  %subplot(1,4,i);
  figure
  for jj=1:length(col_p)
      filtered=[];
      [s,filtered]=medoutlierfilt_nan(hg1(:,col_p(jj)),2,'false',1);
      hg1(:,col_p(jj))=filtered;
      subplot(3,2,jj);
      plot(hg1(:,1),hg1(:,col_p(jj)),'.');
      axis('tight')
      datetick('x',2,'keeplimits','keepticks')
      title(hg_leg(col_p(jj)))
  end
  suptitle(Cal.brw_str(i))
end

%% Ploteos de sl_cr
close all; % !Un poco antes de Strong wind hay cambio (ultima entrada)!
config_events.limits={[330 370],[300 360],[290 330]};
if exist('sl_cr','var')
    for i=1:Cal.n_brw
%         fid=fopen(['log',num2str(Cal.brw(i)),'.txt'],'r'); log=textscan(fid,'%f,%q');

        figure;
        aux=cell2mat(sl_cr{i}); hg_idx=~aux(:,2)==0;
        plot(aux(hg_idx,1),aux(hg_idx,22),'.');
        hold on; plot(aux(~hg_idx,1),aux(~hg_idx,22),'xr');
        set(gca,'YLim',config_events.limits{i}); title(['Standard Lamp R6   ',Cal.brw_name{i}]);
%         vline(log{1}(2:end),'-k',log{2}(2:end))
        datetick('x',26,'Keeplimits','Keepticks'); grid;  %set(gca,'YLim',[300 380]); 
    end
end

%% SL Report
close all;
f0=figure;
sl_s={}; slf={}; R6_={};

for ii=1:3
    try
    if ii==Cal.n_inst 
      [slf{ii},sl_s{ii},sl_out{ii},R6_{ii}]=sl_report_jday(ii,sl,Cal.brw_str,...
                               'date_range',datenum(Cal.Date.cal_year,1,1),...
                               'outlier_flag',1,'fplot',1);
    else
      [slf{ii},sl_s{ii},sl_out{ii},R6_{ii}]=sl_report_jday(ii,sl,Cal.brw_str,...
                               'date_range',datenum(Cal.Date.cal_year,1,1),...
                               'outlier_flag',1,'fplot',1);
    end
    catch
     aux=lasterror;aux.message  
     disp(Cal.brw_str(i)) 
    end
end
% anadimos el sl
save(Cal.file_save,'-APPEND','slf','sl_s','sl_out','R6_')

%% READ Configuration
% configuration for FINAL days
close all
A=[]; ETC=[]; SL_B=[]; cfg=[]; icf_brw=[]; 

[A,ETC,SL_B,cfg,icf_brw]=read_cal_config(config,Cal,sl_s);

try
 tabla_sl=printmatrix([Cal.brw(:)';ETC.old;ETC.new;A.old;A.new;Cal.SL_OLD_REF';Cal.SL_NEW_REF';fix(SL_B)']',4);
 cell2csv(fullfile(Cal.file_latex,'tabla_config.csv'),tabla_sl',';');
catch
  l=lasterror;
  disp(l.message)
  disp('tabla de configuracion posiblemente incompleta')
end

matrix2latex(round(10000*[Cal.brw(:)';ETC.old;ETC.new;A.old;A.new;Cal.SL_OLD_REF';Cal.SL_NEW_REF';fix(SL_B)']')/10000,...
    fullfile(Cal.file_latex,['table_SL_',Cal.brw_str{Cal.n_inst},'.tex']));

%% Data recalculation for summaries  and individual observations
cal=cell(1,Cal.n_brw);
summary=cal;
summary_old=cal;
refs.SL_NEW_REF=Cal.SL_NEW_REF;
refs.SL_OLD_REF=Cal.SL_OLD_REF;

for i=1:Cal.n_brw
    if i==Cal.n_inst
     [cal{i},summary{i},summary_old{i}]=summary_reprocess(refs,i,ozone_ds,ozone_sum,A,sl_s,1);
    else
     [cal{i},summary{i},summary_old{i}]=summary_reprocess(refs,i,ozone_ds,ozone_sum,A,sl_s,1);
    end
end

% para detectar superoutliers
figure; plot(summary{1}(:,1),summary{1}(:,6),'.',...
             summary{2}(:,1),summary{2}(:,6),'.',...
             summary{3}(:,1),summary{3}(:,6),'.');
grid, legend(gca,'ref 157','ref 183','ref 185','Location','SouthEast'); datetick('x',26);

% id_out=find(summary{1}(:,6)>370 | summary{1}(:,6)<220);
% disp('#157'); disp(unique(diaj(summary{1}(id_out,1))));
% summary{1}(id_out,:)=[]; summary_old{1}(id_out,:)=[];
% id_out=find(summary{3}(:,6)>370 | summary{3}(:,6)<220);
% disp('#185'); disp(unique(diaj(summary{3}(id_out,1))));
% summary{3}(id_out,:)=[]; summary_old{3}(id_out,:)=[];
% id_out=find(summary{2}(:,6)>370 | summary{2}(:,6)<220);
% disp('#183'); disp(unique(diaj(summary{2}(id_out,1))));
% summary{2}(id_out,:)=[]; summary_old{2}(id_out,:)=[];
% 
% figure; plot(summary{1}(:,1),summary{1}(:,10),'.',...
%              summary{2}(:,1),summary{2}(:,10),'.',...
%              summary{3}(:,1),summary{3}(:,10),'.');
% grid, legend(gca,'ref 157','ref 183','ref 185','Location','SouthEast'); 
% datetick('x',2,'Keeplimits','Keepticks');

summary_orig=summary;
summary_orig_old=summary_old;
save(Cal.file_save,'-APPEND','summary_old','summary_orig_old','summary','summary_orig');

%% Filter distribution 
for i=1:3
    f=figure;
[aux,b]=histc(summary_old{i}(:,5),[0,64,128,192,256,512]);
label_=mmcellstr(sprintf('Filter #%d=  %.1f%% |',[0:5;100*aux'/sum(aux)]));
h=pie3(aux,label_);
title(sprintf('%s%s','Filter distribution Measurements,  ',Cal.brw_name{i}));
set(f,'Tag','FILTER_DISTRIBUTION');
end

% %% Si queremos eliminar algun filtro
%  j=find(summary{2}(:,5)==256); summary{2}(j,:)=[];
%  j=find(summary_old{2}(:,5)==256); summary_old{2}(j,:)=[];

%% FILTER ANALISYS
%[etc,media,fi,fi_avg]=filter_rep('185','date_range',datenum(2010,09,01));
%
%

%% filter correction 185

ETC_C=[0,0,0,7,15,0];
[summary_old_corr summary_corr]=filter_corr(summary_orig,summary_orig_old,Cal.n_inst,A,ETC_C);
summary_old{Cal.n_inst}=summary_old_corr; summary{Cal.n_inst}=summary_corr;

% %% Depuración de outliers
% [outlier,data_out]=summary_dep(summary_old{Cal.n_inst},'lim',2);
% [outlier,data_out]=summary_dep(inst_,nref_,2,'lim',2);

%% Trabajamos con la configuración del #157 original, icf28406.157,
n_ref=Cal.n_ref(1); % Cuidado con cual referencia estamos asignando
n_inst=Cal.n_inst;
brw=Cal.brw;

brw_str=Cal.brw_str;
blinddays={};

%%
close all
x={};  ox={}; osc_smooth={};
blinddays={};  blinddays_={};  blinddays__={};
 
for i=1:Cal.n_brw
    blinddays{i}=[300:365];
end;

for i=1:Cal.n_brw
    blinddays_{i}=[001:100]; 
end;

jday=findm(diaj(summary{Cal.n_inst}(:,1)),blinddays{Cal.n_inst},0.5);
n185=summary{Cal.n_inst}(jday,:);
jday_=findm(diaj(summary{Cal.n_inst}(:,1)),blinddays_{Cal.n_inst},0.5);
n185_=summary{Cal.n_inst}(jday_,:);

% ref
ref=Cal.n_ref(1); jday_ref=findm((diaj(summary{ref}(:,1))),blinddays{ref},0.5);
n157=summary{ref}(jday_ref,:);
ref=Cal.n_ref(1); jday_ref_=findm((diaj(summary{ref}(:,1))),blinddays_{ref},0.5);
n157_=summary{ref}(jday_ref_,:);

% 183
    jday_183=findm((diaj(summary{2}(:,1))),blinddays{2},0.5);
n183=summary{2}(jday_183,:);
    jday_183_=findm((diaj(summary{2}(:,1))),blinddays_{2},0.5);
n183_=summary{2}(jday_183_,:);

%%
 [r,df,r185_2010,dat]=ratio_min(n185(:,[1,6,3,2,8,9,4,5]),n157(:,[1,6,3,2,8,9,4,5]),2);
 [r,df,r185_2011,dat]=ratio_min(n185_(:,[1,6,3,2,8,9,4,5]),n157_(:,[1,6,3,2,8,9,4,5]),2);
     r_aux=cat(1,r185_2010,r185_2011);
        rat=figure; ha=tight_subplot(3,1,.05,.1);    
        axes(ha(1)); [rp_grp mea_n]=plot_ratio(r_aux,'inst',Cal.brw_str{3},'ref',Cal.brw_str{ref}); 
        axes(ha(1)); set(gca,'YLim',[-2.5 2.5],'XTickLabel',[]);  grid
%         h=vline(datejuli(2010,265),'-r'); set(h,'LineWidth',2); 

 [r,df,r183_2010,dat]=ratio_min(n183(:,[1,6,3,2,8,9,4,5]),n157(:,[1,6,3,2,8,9,4,5]),2);
 [r,df,r183_2011,dat]=ratio_min(n183_(:,[1,6,3,2,8,9,4,5]),n157_(:,[1,6,3,2,8,9,4,5]),2);
     r_aux=cat(1,r183_2010,r183_2011);
        axes(ha(2)); [rp_grp mea_n]=plot_ratio(r_aux,'inst',Cal.brw_str{2},'ref',Cal.brw_str{ref}); 
        axes(ha(2)); set(gca,'YLim',[-2.5 2.5],'XTickLabel',[]);  grid
%         h=vline(datejuli(2010,265),'-r'); set(h,'LineWidth',2); 
        legend('off');
       
 [r,df,rcros_2010,dat]=ratio_min(n183(:,[1,6,3,2,8,9,4,5]),n185(:,[1,6,3,2,8,9,4,5]),2);
 [r,df,rcros_2011,dat]=ratio_min(n183_(:,[1,6,3,2,8,9,4,5]),n185_(:,[1,6,3,2,8,9,4,5]),2);
     r_aux=cat(1,rcros_2010,rcros_2011); 
        axes(ha(3)); plot_ratio(r_aux,'inst',Cal.brw_str{2},'ref',Cal.brw_str{3});
        axes(ha(3)); set(gca,'YLim',[-2.5 2.5]);  legend('off'); grid
%         h=vline(datejuli(2010,265),'-r'); set(h,'LineWidth',2); 
        datetick('x',2,'Keeplimits','Keepticks'); 


%%        
%  [r,df,r185_2009orig,dat]=ratio_min(inst0(:,[1,6,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),2);
     r_aux=cat(1,r185_2010,r185_2011);
        rat=figure; ha=tight_subplot(2,1,.05,.1);    
        axes(ha(1)); [rp_grp mea_n]=plot_ratio(r_aux,'inst',brw_str{3},'ref',brw_str{ref}); 
        axes(ha(1)); set(gca,'YLim',[-2.5 2.5],'XTickLabel',[]);  grid
        h=vline([datejuli(2009,245),datejuli(2009,265)],'-k'); set(h,'LineWidth',2);
        vline(datejuli(2009,204),'k'); text(r_aux(1,1)-20,2,'\bfa)');

figure; set(gcf,'tag','RATIO_NEWs'); n_w=tight_subplot(2,1,.05,.1);
axes(n_w(1)); plot(r_aux(:,1),r_aux(:,2),'*'); 
hold on;
g=gscatter(rp_grp(:,1),rp_grp(:,2),rp_grp(:,3));
legend(g,'location','North','orientation','horizontal');
set(g,'Marker','.','MarkerSize',5); title([brw_str{Cal.n_inst},' - ',brw_str{Cal.n_ref(1)},' / ',brw_str{Cal.n_ref(1)},  '    RBCCE config']);
set(gca,'Ylim',[-3.5 3.5],'XTickLabel',[]);  grid    

[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary{1},summary{2},summary{3}},0.3410,...
                       Cal.n_inst,Cal.n_ref(1),2,1.6,0.03,blinddays_{3});
% o3r= (n185_(:,8)-1580)./(0.3410*n185_(:,3)*10);%blinddays{3}
% n185_(:,12)=o3r;                
%     [x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
%        n185_(:,[1,10,3,2,8,9,4,5]),n157_(:,[1,6,3,2,8,9,4,5]),...
%        2,Cal.brw_str{3},Cal.brw_str{ref});
[r,df,r185_2010rbcc,dat]=ratio_min(inst_(:,[1,12,3,2,8,9,4,5]),nref_(:,[1,6,3,2,8,9,4,5]),2);

% El 265 hay cambio
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary{1},summary{2},summary{3}},0.3410,...
                       Cal.n_inst,Cal.n_ref(1),2,1.6,0.03,blinddays__{3});
o3r= (inst__(:,8)-1580)./(0.3410*inst__(:,3)*10);
inst__(:,12)=o3r;                
    [x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       inst__(:,[1,6,3,2,8,9,4,5]),nref__(:,[1,6,3,2,8,9,4,5]),...
       2,Cal.brw_str{3},Cal.brw_str{ref});
[r,df,r185_2010rbcc2,dat]=ratio_min(inst__(:,[1,12,3,2,8,9,4,5]),nref__(:,[1,6,3,2,8,9,4,5]),2);

% Llevamos resultados a ASCII
aux=cat(1,r185_2010rbcc,r185_2010rbcc2);
fid=fopen('ESA2010.txt','w');
fprintf(fid,'%s   %s    %s   %s    %s   %s    %s   %s\r\n',...
    'date(ref)', 'rel.dif.', 'sza(ref)', 'm(ref)', 'ozono(inst)',...
    'ozono(ref)', 'temp(inst)', 'filter(inst)');
for i=1:size(aux,1)
    fprintf(fid,'%f   %f    %f   %f    %f   %f    %f   %f\r\n',aux(i,1),aux(i,2),...
            aux(i,3),aux(i,4),aux(i,5),aux(i,6),aux(i,7),aux(i,8));
end
fclose all;
A=textscan(fopen('ESA2010.txt','r'),'%f %f %f %f %f %f %f %f','Headerlines',1);
inst_(find(inst_(:,12)<0),:)=[];
figure; plot(A{1},A{6},'r*');
hold on; plot(A{1},A{5},'+');
hold on; plot(cat(1,inst_(:,1),inst__(:,1)),cat(1,inst_(:,12),inst__(:,12)),'g.');
hold on; plot(cat(1,nref_(:,1),nref__(:,1)),cat(1,nref_(:,6),nref__(:,6)),'b.');
legend(gca,'sync inst','sync ref','all inst','all ref');
datetick('x',2,'Keeplimits','Keepticks');

fid=fopen('Brw185_ESA2010.txt','w');
fprintf(fid,'%s   %s    %s   %s    %s   %s    %s  %s\r\n',...
    'date','sza', 'm','temp','filt','ozone', 'std','MS9');
n185=cat(1,inst_,inst__);
for i=1:size(n185,1)
  fprintf(fid,'%f %f %f %f %f %f %f %f\r\n',n185(i,1),n185(i,2),n185(i,3),n185(i,4),...
                                         n185(i,5),n185(i,12),n185(i,13),n185(i,8));
end
fclose all;

fid=fopen('Brw157_ESA2010.txt','w');
fprintf(fid,'%s   %s    %s   %s    %s   %s    %s    %s\r\n',...
    'date','sza', 'm','temp','filt','ozone', 'std', 'MS9');
n157=cat(1,nref_,nref__);
for i=1:size(n157,1)
  fprintf(fid,'%f %f %f %f %f %f %f %f\r\n',n157(i,1),n157(i,2),n157(i,3),n157(i,4),...
                                            n157(i,5),n157(i,6),n157(i,7),n157(i,8));
end
fclose all;



%%
% A partir del 20409 hay un cambio (ver arriba). 
% Pero para calibrar uso desde el 15/Agosto (22709), para no trabajar con datos anomalos
% (no hay cambio debido a Arenosillo2009)
figure;
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary{1},summary{2},summary{3}},0.3415,...% A.new(Cal.n_inst)
                       Cal.n_inst,Cal.n_ref(1),2,2,0.03,[227:245 265:360]);% [227:245 265:296 298:351 354:360]   
o3r= (inst_(:,8)-1595)./(0.3415*inst_(:,3)*10);
inst_(:,12)=o3r;                
    [x,r185_2009rbcc_,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       inst_(:,[1,12,3,2,8,9,4,5]),nref_(:,[1,6,3,2,8,9,4,5]),...
       2,brw_str{3},brw_str{ref});
[r,df,r185_2009rbcc_,dat]=ratio_min(inst_(:,[1,12,3,2,8,9,4,5]),nref_(:,[1,6,3,2,8,9,4,5]),2);

     r_aux=cat(1,r185_2009orig,r185_2009rbcc,r185_2009rbcc2,r185_2009rbcc_);
        axes(ha(2)); [rp_grp mea_n]=plot_ratio(r_aux,'inst',brw_str{3},'ref',brw_str{ref}); 
        axes(ha(2)); set(gca,'YLim',[-2.5 2.5],'XTickLabel',[]);  grid
        h=vline([datejuli(2009,245),datejuli(2009,265)],'-k'); set(h,'LineWidth',2);
        h=vline(datejuli(2009,204),'-r','zn #157'); set(h,'LineWidth',2);
        h=vline(datejuli(2009,85),'-r','Storm power off'); set(h,'LineWidth',2);
        datetick('x',2,'Keeplimits','Keepticks'); legend('off'); grid
                
        axes(n_w(2)); plot(r_aux(:,1),r_aux(:,2),'*'); 
        hold on;
        g=gscatter(rp_grp(:,1),rp_grp(:,2),rp_grp(:,3));
        set(g,'Marker','.','MarkerSize',5);    set(gca,'Ylim',[-3.5 3.5]);    
        legend('off'); title(''); grid
        datetick('x',2,'Keeplimits','Keepticks'); 

% disp('mean month ratio statistics');
% tableform({'orig', 'orig_std', 'final', 'final_std'},...
%    [mea_n.media(:,5), mea_n.sigma(:,5), mea_new.media(:,5), mea_new.sigma(:,5)]);

% Llevamos resultados a ASCII
fid=fopen('are09.txt','w');
fprintf(fid,'%s   %s    %s   %s    %s   %s    %s   %s\r\n',...
    'date(ref)', 'rel.dif.', 'sza(ref)', 'm(ref)', 'ozono(inst)',...
    'ozono(ref)', 'temp(inst)', 'filter(inst)');
for i=1:size(r185_2009rbcc_,1)
    fprintf(fid,'%f   %f    %f   %f    %f   %f    %f   %f\r\n',r185_2009rbcc_(i,1),r185_2009rbcc_(i,2),...
            r185_2009rbcc_(i,3),r185_2009rbcc_(i,4),r185_2009rbcc_(i,5),r185_2009rbcc_(i,6),r185_2009rbcc_(i,7),r185_2009rbcc_(i,8));
end
fclose all;
A=textscan(fopen('are09.txt','r'),'%f %f %f %f %f %f %f %f','Headerlines',1);
inst_(find(inst_(:,12)<0),:)=[];
figure; plot(A{1},A{6},'r*')
hold on; plot(A{1},A{5},'+')
hold on; plot(inst_(:,1),inst_(:,12),'g.')
hold on; plot(nref_(:,1),nref_(:,6),'b.')


fid=fopen('Brw185_are09.txt','w');
fprintf(fid,'%s   %s    %s   %s    %s   %s    %s\r\n',...
    'date','sza', 'm','temp','filt','ozone', 'std');
for i=1:size(inst_,1)
  fprintf(fid,'%f %f %f %f %f %f %f\r\n',inst_(i,1),inst_(i,2),inst_(i,3),inst_(i,4),inst_(i,5),inst_(i,12),inst_(i,13));
end
fclose all;

fid=fopen('Brw157_are09.txt','w');
fprintf(fid,'%s   %s    %s   %s    %s   %s    %s\r\n',...
    'date','sza', 'm','temp','filt','ozone', 'std');
for i=1:size(nref_,1)
  fprintf(fid,'%f %f %f %f %f %f %f\r\n',nref_(i,1),nref_(i,2),nref_(i,3),nref_(i,4),nref_(i,5),nref_(i,6),nref_(i,7));
end
fclose all;


%% 183: Punto de partida: ratios con las calibraciónes originales (icf18008 hasta IOS09, a partir de ahi icf24709).

    r_aux=cat(1,r183_2009orig,r183_2009orig_1,r183_2009orig_,r183_2009orig_ios);% ,r185_2009ios__,r185_2009ios_ 
        rat=figure; ha=tight_subplot(2,1,.05,.1);    
        axes(ha(1)); [rp_grp mea_n]=plot_ratio(r_aux,'inst',brw_str{3},'ref',brw_str{ref}); 
        axes(ha(1)); set(gca,'YLim',[-2.5 2.5],'XTickLabel',[]);  grid
        h=vline([datejuli(2009,245),datejuli(2009,265)],'-k'); set(h,'LineWidth',2);
        h=vline(datejuli(2009,204),'-r'); set(h,'LineWidth',2); 
        h=vline(datejuli(2009,85),'-r','Storm power off'); set(h,'LineWidth',2);
 
figure; set(gcf,'tag','RATIO_NEWs'); n_w=tight_subplot(2,1,.05,.1);
axes(n_w(1)); plot(r_aux(:,1),r_aux(:,2),'*'); 
hold on;
g=gscatter(rp_grp(:,1),rp_grp(:,2),rp_grp(:,3));
legend(g,'location','North','orientation','horizontal');
set(g,'Marker','.','MarkerSize',5); title([brw_str{Cal.n_inst},' - ',brw_str{Cal.n_ref(1)},' / ',brw_str{Cal.n_ref(1)},  '    RBCCE config']);
set(gca,'Ylim',[-3.5 3.5],'XTickLabel',[]);  grid    

% Recalibramos el 183 en el periodo [1:078]. El 079 cambio gel silice y
% problemas con zn
% No puedo mejorar lo que hay. Revisar Temp. Coeffs. TRabajo con originales (0)
figure;
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary_old{1},summary_old{2},summary_old{3}},0.3415,...% A.new(Cal.n_inst)
                       2,Cal.n_ref(1),1,1.6,0.03,1:78);%   
jday_183_=findm((diaj(summary_old{2}(:,1))),1:78,0.5);
old_n183=summary_old{2}(jday_183_,:);

o3r= (old_n183(:,8)-1725)./(0.3415*old_n183(:,3)*10);
old_n183(:,12)=o3r;                
    [x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       old_n183(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),...
       2,brw_str{2},brw_str{ref});
[r,df,r183_raro_1,dat]=ratio_min(old_n183(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),2);

% Mantenemos icf18008 en el periodo [100:195] (recalibrando coincide, pero llamamos icf08009.183)
figure;
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary_old{1},summary_old{2},summary_old{3}},0.3415,...% A.new(Cal.n_inst)
                       2,Cal.n_ref(1),1,1.6,0.03,[101:120 121 122 132:144 146 147 153:187 189:195]);%   
jday_183_=findm((diaj(summary_old{2}(:,1))),[080:120 121 122 132:144 146 147 153:187 189:195],0.5);
old_n183=summary_old{2}(jday_183_,:);

o3r= (old_n183(:,8)-1735)./(0.3415*old_n183(:,3)*10);
old_n183(:,12)=o3r;                
    [x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       old_n183(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),...
       2,brw_str{2},brw_str{ref});
[r,df,r183_raro_2,dat]=ratio_min(old_n183(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),2);
   

% Existe Cambio apreciable en la respuesta del # desde el día 195 hasta IOS09
% No se puede recalibrar (periodo [196:247]). ¿Motivos? Nada apuntado en log
% La respuesta del equipo cambia rapidamente en 1 mes aproximadamente (amplitud ~3%)
% De abajo se distinguen tres periodos: [198:210], [220:233] y [234:246]
    [x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       n183_(:,[1,10,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),...
       2,brw_str{2},brw_str{ref});

figure;
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary_old{1},summary_old{2},summary_old{3}},0.3415,...% A.new(Cal.n_inst)
                       2,Cal.n_ref(1),1,1,0.03,200:210);%   
jday_183_=findm((diaj(summary_old{2}(:,1))),196:210,0.5);
old_n183_=summary_old{2}(jday_183_,:);

o3r= (old_n183_(:,8)-1745)./(0.3415*old_n183_(:,3)*10);
old_n183_(:,12)=o3r;                
    [x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       old_n183_(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),...
       2,brw_str{2},brw_str{ref});
[r,df,r183__1,dat]=ratio_min(old_n183_(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),2);

figure;
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary_old{1},summary_old{2},summary{3}},0.3415,...% A.new(Cal.n_inst)
                       2,Cal.n_ref(1),1,1,0.03,220:233);%   
jday_183_=findm((diaj(summary_old{2}(:,1))),210:233,0.5);
old_n183_=summary_old{2}(jday_183_,:);

o3r= (old_n183_(:,8)-1755)./(0.3415*old_n183_(:,3)*10);
old_n183_(:,12)=o3r;                
    [x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       old_n183_(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),...
       2,brw_str{2},brw_str{ref});
[r,df,r183__2,dat]=ratio_min(old_n183_(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),2);

figure;
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary_old{1},summary_old{2},summary{3}},0.3415,...% A.new(Cal.n_inst)
                       2,Cal.n_ref(1),1,1,0.03,234:245);%   
jday_183_=findm((diaj(summary_old{2}(:,1))),234:245,0.5);
old_n183_=summary_old{2}(jday_183_,:);

o3r= (old_n183_(:,8)-1770)./(0.3415*old_n183_(:,3)*10);
old_n183_(:,12)=o3r;                
    [x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       old_n183_(:,[1,12,3,2,8,9,4,5]),nref_(:,[1,10,3,2,8,9,4,5]),...
       2,brw_str{2},brw_str{ref});
[r,df,r183__3,dat]=ratio_min(old_n183_(:,[1,12,3,2,8,9,4,5]),nref0(:,[1,6,3,2,8,9,4,5]),2);
                
% Pasamos a evaluar actuación de IOS09, hasta final de año
figure;
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary{1},summary_old{2},summary{3}},0.3415,...% A.new(Cal.n_inst)
                       2,Cal.n_ref(1),1,1.4,0.03,248:264);%   
jday_183_=findm((diaj(summary_old{2}(:,1))),248:264,0.5);
n183_parc=summary_old{2}(jday_183_,:);
o3r= (n183_parc(:,8)-1625)./(0.3415*n183_parc(:,3)*10);
n183_parc(:,12)= o3r;   
[x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       n183_parc(:,[1,12,3,2,8,9,4,5]),nref_(:,[1,6,3,2,8,9,4,5]),...
       2,brw_str{2},brw_str{ref});               
[r,df,r183_parc,dat]=ratio_min(n183_parc(:,[1,12,3,2,8,9,4,5]),nref_(:,[1,6,3,2,8,9,4,5]),2);

figure;
[ETC_BLIND,o3c_BLIND,m_etc_BLIND]=ETC_calibration_C(Cal,{summary{1},summary_old{2},summary{3}},0.3415,...% A.new(Cal.n_inst)
                       2,Cal.n_ref(1),1,1.4,0.03,blinddays_ios{2});%   
o3r= (n183_ios(:,8)-1630)./(0.3415*n183_ios(:,3)*10);
n183_ios(:,12)= o3r;   
[x,r,rp,ra,dat,ox,osc_smooth185_rbcc_]=ratio_min_ozone(...
       n183_ios(:,[1,12,3,2,8,9,4,5]),nref_(:,[1,6,3,2,8,9,4,5]),...
       2,brw_str{2},brw_str{ref});               
[r,df,r183_2009_ios,dat]=ratio_min(n183_ios(:,[1,12,3,2,8,9,4,5]),nref_(:,[1,6,3,2,8,9,4,5]),2);

  r_aux=cat(1,r183_2009orig,r183_2009orig_1,r183_2009orig_,r183_2009orig_ios,...
         r183_raro_1,r183_raro_2,r183__1,r183__2,r183__3,r183_2009_ios);% ,r185_2009ios__,r185_2009ios_ 
        axes(ha(2)); [rp_grp mea_n]=plot_ratio(r_aux,'inst',brw_str{3},'ref',brw_str{ref}); 
        axes(ha(2)); set(gca,'YLim',[-2.5 2.5]); title(''); grid
        h=vline([datejuli(2009,245),datejuli(2009,265)],'-k'); set(h,'LineWidth',2);
        h=vline(datejuli(2009,204),'-r','zn #157'); set(h,'LineWidth',2);
        h=vline(datejuli(2009,85),'-r','Storm power off'); set(h,'LineWidth',2);
        datetick('x',2,'Keeplimits','Keepticks'); legend('off'); grid

axes(n_w(2)); plot(r_aux(:,1),r_aux(:,2),'*'); 
hold on;
g=gscatter(rp_grp(:,1),rp_grp(:,2),rp_grp(:,3));
legend(g,'location','North','orientation','horizontal');
set(g,'Marker','.','MarkerSize',5); title([brw_str{2},' - ',brw_str{Cal.n_ref(1)},' / ',brw_str{Cal.n_ref(1)},  '    RBCCE config']);
set(gca,'Ylim',[-5 5],'XTickLabel',[]);  legend('off'); title(''); grid    
datetick('x',2,'Keeplimits','Keepticks'); 

%% Ratios cruzadas
% Resumen:
figure
 [r,df,rcros_2009orig_1,dat]=ratio_min(inst(:,[1,10,3,2,8,9,4,5]),n183_old(:,[1,6,3,2,8,9,4,5]),2);
 [r,df,rcros_2009orig_2,dat]=ratio_min(inst(:,[1,10,3,2,8,9,4,5]),n183_1(:,[1,10,3,2,8,9,4,5]),2);
 [r,df,rcros_2009orig_3,dat]=ratio_min(inst_(:,[1,10,3,2,8,9,4,5]),n183_(:,[1,10,3,2,8,9,4,5]),2);
 [r,df,rcros_2009orig_4,dat]=ratio_min(inst_(:,[1,10,3,2,8,9,4,5]),n183_ios(:,[1,6,3,2,8,9,4,5]),2);

     r_aux=cat(1,rcros_2009orig_1,rcros_2009orig_2,rcros_2009orig_3,rcros_2009orig_4); % ,rcros_2009ios__
        figure; ha=tight_subplot(2,1,.05);
        axes(ha(1)); plot_ratio(r_aux,'inst',brw_str{2},'ref',brw_str{3});
        set(gca,'YLim',[-2.5 2.5],'XtickLabel',[]); grid
 
 [r,df,rcros_2009_1,dat]=ratio_min(inst(:,[1,12,3,2,8,9,4,5]),n183_old(:,[1,6,3,2,8,9,4,5]),2);
 [r,df,rcros_2009_2,dat]=ratio_min(inst(:,[1,12,3,2,8,9,4,5]),n183_1(:,[1,12,3,2,8,9,4,5]),2);
 [r,df,rcros_2009_3,dat]=ratio_min(inst_(:,[1,10,3,2,8,9,4,5]),n183_(:,[1,10,3,2,8,9,4,5]),2);
 [r,df,rcros_2009_4,dat]=ratio_min(inst_(:,[1,12,3,2,8,9,4,5]),n183_ios(:,[1,12,3,2,8,9,4,5]),2);

     r_aux=cat(1,rcros_2009_1,rcros_2009_2,rcros_2009_3,rcros_2009_4); % ,rcros_2009ios__
        axes(ha(2)); 
        plot_ratio(r_aux,'inst',brw_str{2},'ref',brw_str{3});
        set(gca,'YLim',[-2.5 2.5]);  legend('off'); title(''); grid
        datetick('x',2,'Keeplimits','Keepticks'); 


%% para chequear matrices de configuración
