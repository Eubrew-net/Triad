clear all    
addpath(genpath(fullfile('~','CODE','rbcce.aemet.es','iberonesia','matlab')));
path_root=(fullfile('~','CODE','rbcce.aemet.es','iberonesia','RBCC_E','Triad'))


%Cargamos los datos/a�os y los agrupamos por brewer en la misma variable.
             
aux=cell(3,1);
for ano=2019:-1:2005
    for jj=1:3
       brw=[157,183,185]; 
       filename=fullfile(path_root,num2str(ano),'Triad',sprintf('summary_old_Brw%03d_%04d.txt',brw(jj),ano))
       summary_old=load(filename);
       inicio=datenum(ano,1,1);fin=datenum(ano,1,366);
       b=find(summary_old(:,1)>=inicio & summary_old(:,1)<=fin);
       c=summary_old(b,:);
       b=[c(:,1),c(:,2),c(:,3),c(:,6),c(:,12)]; 
       aux{jj}=[b;aux{jj}];
    end
end

aux157oma=aux{1}; aux183oma=aux{2}; aux185oma=aux{3};
aux157=aux{1}; aux183=aux{2}; aux185=aux{3};

%% Datos para Omaira y comparativa entre ALberto, eubrewnet y mis datos
% 
% % En esta secci�n preparo los datos que le voy a dar a Omaira. Para ello,
% % al comienzo del codigo guarde las variables auxXXXoma. 
% 
 B{1}=aux157oma; B{2}=aux183oma; B{3}=aux185oma;
 BXXXoma{1}=[]; BXXXoma{2}=[]; BXXXoma{3}=[]; 
 
 for j=1:1:3
     a=size(B{j});
     BXXXoma{j}=nan(a(1,1),8);
     
     fecha=datevec(B{j}(:,1)); BXXXoma{j}(:,1)=fecha(:,1); % a�o
     fecha0=fecha; fecha0(:,2)=1; fecha0(:,3)=1;
     fecha0=fix(datenum(fecha0)); BXXXoma{j}(:,2)=fix(B{j}(:,1)-fecha0+1); % D�a
     BXXXoma{j}(:,3)=fecha(:,4)*60+fecha(:,5);  % Hora medida (minutos)
     BXXXoma{j}(:,4)=B{j}(:,4); %ozono 
     BXXXoma{j}(:,5)=B{j}(:,5); %ozono_slcorrec 
     BXXXoma{j}(:,6)=B{j}(:,2); %SZA
     BXXXoma{j}(:,7)=B{j}(:,3); %masa �ptica
     BXXXoma{j}(:,8)=0; % Variable de control
 end
 
 aux157oma=BXXXoma{1}; aux183oma=BXXXoma{2}; aux185oma=BXXXoma{3};
% % FILTRADO: Eliminamos las medidas realizadas en campa�as de calibraci�n 
% % o cuanto el equipo esta averiado. 
%%Brewer#185
%
aa=[];
a=find(aux185oma(:,1)==2007 & aux185oma(:,2)>=245 & aux185oma(:,2)<=254); aa=[aa;a];
a=find(aux185oma(:,1)==2009 & aux185oma(:,2)>=246 & aux185oma(:,2)<=265); aa=[aa;a];
a=find(aux185oma(:,1)==2010 & aux185oma(:,2)>=200 & aux185oma(:,2)<=212); aa=[aa;a];
a=find(aux185oma(:,1)==2011 & aux185oma(:,2)>=186 & aux185oma(:,2)<=197); aa=[aa;a];
a=find(aux185oma(:,1)==2012 & aux185oma(:,2)>=164 & aux185oma(:,2)<=181); aa=[aa;a];
a=find(aux185oma(:,1)==2013 & aux185oma(:,2)>=155 & aux185oma(:,2)<=176); aa=[aa;a];
a=find(aux185oma(:,1)==2014 & aux185oma(:,2)>=184 & aux185oma(:,2)<=209); aa=[aa;a];
a=find(aux185oma(:,1)==2016 & aux185oma(:,2)>=142 & aux185oma(:,2)<=157); aa=[aa;a];      
a=find(aux185oma(:,1)==2005 & aux185oma(:,2)==276); aa=[aa;a];
a=find(aux185oma(:,1)==2006 & aux185oma(:,2)==100); aa=[aa;a];
a=find(aux185oma(:,1)==2006 & aux185oma(:,2)==30);  aa=[aa;a];
a=find(aux185oma(:,1)==2006 & aux185oma(:,2)==323); aa=[aa;a];
a=find(aux185oma(:,1)==2006 & aux185oma(:,2)==75);  aa=[aa;a];
aux185oma(aa,:)=[];
B{1}(aa,:)=[];
%%Brewer#183
aa=[];
a=find(aux183oma(:,1)==2013 & aux183oma(:,2)>=155 & aux183oma(:,2)<=176); aa=[aa;a];       
a=find(aux183oma(:,1)==2006 & aux183oma(:,2)>=180& aux183oma(:,2)<=282); aa=[aa;a]; 
aux183oma(aa,:)=[];
B{2}(aa,:)=[]
%%Brewer#157
aa=[];
a=find(aux157oma(:,1)==2006 & aux157oma(:,2)==160); 
aux157oma(a,:)=[]; B{3}(a,:)=[]
a=find(aux157oma(:,1)==2014 & aux157oma(:,2)>=130 & aux157oma(:,2)<=161); 
aux157oma(a,:)=[]; B{3}(a,:)=[]

% Debemos seleccionar los periodos donde se emplea el ozono o ozono_sl correc
B{1}=[B{1},(ones(size(B{1},1),1)-1)]
B{2}=[B{2},(ones(size(B{2},1),1)-1)]
B{3}=[B{3},(ones(size(B{3},1),1)-1)]

% Brewer 157
a=find(aux157oma(:,1)==2014 & aux157oma(:,2)>=258 & aux157oma(:,2)<=366);
aux157oma(a,8)=1; B{1}(a,6)=1;
a=find(aux157oma(:,1)==2016 & aux157oma(:,2)>=1 & aux157oma(:,2)<=99);
aux157oma(a,8)=1; B{1}(a,6)=1;

% Brewer 183 
a=find(aux183oma(:,1)==2006 & aux183oma(:,2)>=1 & aux183oma(:,2)<=366);  aux183oma(a,8)=1; B{2}(a,6)=1;
a=find(aux183oma(:,1)==2007 & aux183oma(:,2)>=179 & aux183oma(:,2)<=366); aux183oma(a,8)=1; B{2}(a,6)=1; 
a=find(aux183oma(:,1)==2008 & aux183oma(:,2)>=1 & aux183oma(:,2)<=366); aux183oma(a,8)=1; B{2}(a,6)=1;
a=find(aux183oma(:,1)==2009 & aux183oma(:,2)>=1 & aux183oma(:,2)<=232); aux183oma(a,8)=1; B{2}(a,6)=1;
a=find(aux183oma(:,1)==2012 & aux183oma(:,2)>=345 & aux183oma(:,2)<=366); aux183oma(a,8)=1; B{2}(a,6)=1;
a=find(aux183oma(:,1)==2013 & aux183oma(:,2)>=1 & aux183oma(:,2)<=366); aux183oma(a,8)=1; B{2}(a,6)=1;
a=find(aux183oma(:,1)==2014 & aux183oma(:,2)>=1 & aux183oma(:,2)<=335); aux183oma(a,8)=1; B{2}(a,6)=1;
    
% Brewer 185
a=find(aux185oma(:,1)==2011 & aux185oma(:,2)>=312 & aux185oma(:,2)<=366); aux185oma(a,8)=1; B{3}(a,6)=1;
a=find(aux185oma(:,1)==2012 & aux185oma(:,2)>=1 & aux185oma(:,2)<=61); aux185oma(a,8)=1; B{3}(a,6)=1;
a=find(aux185oma(:,1)==2014 & aux185oma(:,2)>=44 & aux185oma(:,2)<=343); aux185oma(a,8)=1; B{3}(a,6)=1;

%  Guardamos los datos para Omaira. 
save('ozono157sergio.mat','aux157oma');
save('ozono183sergio.mat','aux183oma');
save('ozono185sergio.mat','aux185oma');

% %Bentor
% Eubrewnet=[]
% for i=2005:1:2016
%     i=num2str(i)
%     a=strcat('C:\CODE\iberonesia\RBCC_Ea\2014\Triad\Langley\OzonosinSL\157_',i,'_ozone_product_1_5.txt');
%     a=load(a);
%     Eubrewnet=[Eubrewnet;a];
% end
% 
% datos=Eubrewnet; Eubrewnet=[];
% Eubrewnet=[datos(:,5),datos(:,6),datos(:,7),datos(:,10),(ones(size(datos,1),1)-1)];
% 
% % Comparamos mis datos y los de Eubrewnet para el 157. Suponiendo que 
% % NUNCA hay que hacer SL correction. (S�lo Brewer 157).
% 
% Comparativa=[]
% for i=1:1:size(B{1},1)
%     t_inicio=B{1}(i,1)-0.00069444444*0.45;
%     t_fin=B{1}(i,1)+0.00069444444*0.45;
%     b=find(Eubrewnet(:,1)>=t_inicio & Eubrewnet(:,1)<=t_fin);
%    if isempty(b)==0
%       c=[B{1}(i,1),B{1}(i,4),Eubrewnet(b(1,1),4)];
%       Comparativa=[Comparativa;c];       
%     end
% end
% 
% ratioSE=Comparativa(:,3)./Comparativa(:,2)
% plot(Comparativa(:,1),ratioSE)
% datetick
% 
% % Comparamos mis datos, Eubrewnet y Alberto ara el 157. Suponiendo que 
% % teniendo en cuenta la SL correction. 
% 
% % Selecciono mis valores de ozono
% 
% OzonoSL=nan(size(B{1},1),4);
% OzonoSL(:,1:3)=B{1}(:,1:3)
% for i=1:1:size(B{1},1)
%     if B{1}(i,6)==0
%        OzonoSL(i,4)=B{1}(i,4);
%     else
%        OzonoSL(i,4)=B{1}(i,5);
%     end
% end
% 
% % Datos de Alberto
% [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13]=textread('O3S0515_ALLSZA.prn',...
%     '%02d:%02d:%02d %d %d %04d %f %f %d %2c %d %f %f ');
% date=datenum(a6,a4,a5,a1,a2,a3);
% jds=strmatch('ds',a10);ty(jds)=1;
% jzs=strmatch('zs',a10);ty(jzs)=2;
% ozo_ds=[date(jds),a13(jds),a7(jds),a12(jds)];
% 
% % Bentor
% 
% EubrewnetSL=[]
% for i=2005:1:2016
%     i=num2str(i)
%     a=strcat('C:\CODE\iberonesia\RBCC_Ea\2014\Triad\Langley\OzonoconSL\157_',i,'_ozone_product_1_5.txt');
%     a=load(a);
%     EubrewnetSL=[EubrewnetSL;a];
% end
% 
% datos=EubrewnetSL; EubrewnetSL=[];
% EubrewnetSL=[datos(:,5),datos(:,6),datos(:,7),datos(:,10),(ones(size(datos,1),1)-1)];
% 
% ComparativaSL=[]
% for i=1:1:size(OzonoSL,1)
%     t_inicio=OzonoSL(i,1)-0.00069444444*0.55;
%     t_fin=OzonoSL(i,1)+0.00069444444*0.55;
%     a=find(ozo_ds(:,1)>=t_inicio & ozo_ds(:,1)<=t_fin);
%     if isempty(a)==0
%         b=find(EubrewnetSL(:,1)>=t_inicio & EubrewnetSL(:,1)<=t_fin);
%         if isempty(b)==0
%             c=[OzonoSL(i,1),OzonoSL(i,4),EubrewnetSL(b(1,1),4),ozo_ds(a(1,1),4)];
%             ComparativaSL=[ComparativaSL;c];
%         end        
%     end
% end
% 
% ratioSESL=ComparativaSL(:,3)./ComparativaSL(:,2)
% ratioSASL=ComparativaSL(:,4)./ComparativaSL(:,2)
% ratioAESL=ComparativaSL(:,3)./ComparativaSL(:,4)
% 
% plot(ComparativaSL(:,1),ratioSESL)
% plot(ComparativaSL(:,1),ratioSASL)
% plot(ComparativaSL(:,1),ratioAESL)
% datetick
% 
% a=importdata('prueba.txt')
% fileID=fopen('prueba.txt')
% 
% C=textread(fileID,...
% '%f %f %f %f %f %f %f %f %f %f %f %f%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f[^\n\r]','Headerlines',52)
% 
% C=textscan(fileID,'%f %f %f','Headerlines',52,'Delimiter',',')
% 
% fid=fopen('prueba.txt')
% fid=fopen(file);
% if fid
%   for i=1:19          % message  
%    c{i}=fgets(fid);
%   end 
%  head=char(c{:});
%  C=textscan(fileID,'%f %16c %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f[^\n\r]','Headerlines',52,'Delimiter',',')
%  fgets(fid);
%  s=fgets(fid)
%  %s=fgets(fid); 
%  %s=fgets(fid);
%  %s=fgets(fid);
%  i=0;
%  config_info={};  % config
%  %nx=8;
%  %nr=4;
%  %while nx<length(s) & nr==4;
%  %i=i+1    
%  try
%  cfg_id=sscanf(s(8:end),'%04d-%02d-%02d, id = %d;');
%  cfg_id=reshape(cfg_id,4,[])
%  catch
%      disp('config not read');
%  end
%  %nx
%  %end
% %  while s(3)=='D'
% %      i=i+1;
% %      config_info{i}=s;
% %      [si,sj]=regexp(s,'[=]\d*');
% %      cfg_id(i)=sscanf(s(si:sj),'=%d');
% %      s=fgets(fid);   % config
% %  end
%  fgets(fid);    
%  for i=1:27 leg{i}=fgets(fid); end%
% 
%  
% fgets(fid);
% format=fgets(fid);
% format=strrep(format,'d','f');  
% format(1)=[];
% head=fgets(fid);
% 
% data=textscan(fid,format);
% data=cell2mat(data);
% 
% %dates to matlab
% fecha_1=datenum(data(:,2:7));
% fecha_p=datenum(data(:,end-5:end));
% fecha_c=datenum(data(:,end-9:end-7));
% 
% data(:,2)=fecha_1;
% data(:,3:7)=[];
% data(:,end-9)=fecha_c;
% data(:,end-8:end-7)=[];
% data(:,end-5)=fecha_p;
% data(:,end-4:end)=[];
% 
% %h=strrep(head,'gmt','YYYY,MM,DD,hh,mm,ss');
% %h=strrep(h,'process_date','pYYYY,pMM,pDD,phh,pmm,pss');
% %h=strrep(h,'configdate','cYYYY,cMM,cDD');
% 
% %date_mat=data(:,10);
% %data=[date_mat,data];
% %data=textscan(fid,'','delimiter','TZ,');
% %date_mat=data(:,7);
% %%salida igual que read_overpas
% %data=[date_mat,data];
% %date_fields=find(~cellfun(@isempty,strfind(leg,'ISO 8601')));
% %date_fileds=[2,27];
% %leg=insertrows(leg',leg(date_fileds)',date_fileds)';
% %date_fileds=[3,28];
% %leg=insertrows(leg',leg(date_fileds)',date_fileds)'; %inser a nan 
% %leg=['Date Matlab',leg(1:3),'NAN',leg(4:end)];  % checl the nan
% %leg=['Date Matlab',leg]; 
% else
%     disp('file error');
% end
% fclose(fid);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %Alberto
% 
% [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13]=textread('O3S0515_ALLSZA.prn',...
%     '%02d:%02d:%02d %d %d %04d %f %f %d %2c %d %f %f ');
% 
% a=size(datos); Eubrewnet=nan(a(1,1),8);
% 
% fecha=datevec(datos(:,5)); Eubrewnet(:,1)=fecha(:,1); % a�o
% fecha0=fecha; fecha0(:,2)=1; fecha0(:,3)=1;
% fecha0=fix(datenum(fecha0)); Eubrewnet(:,2)=fix(datos(:,5)-fecha0+1); % D�a
% Eubrewnet(:,3)=fecha(:,4)*60+fecha(:,5);  % Hora medida (minutos)
% Eubrewnet(:,4)=datos(:,10); %ozono
% Eubrewnet(:,5)=datos(:,10); %ozono_slcorrec
% Eubrewnet(:,6)=datos(:,6); %SZA
% Eubrewnet(:,7)=datos(:,7); %masa �ptica
% Eubrewnet(:,8)=0; % Variable de control