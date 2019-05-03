%% Programa que reproduce el ajuste llevado a cabo por Fieletov para los
% datos medidos en Iza�a por los equipos Brewer.

    % Daily model for midday  
    % O3= A + B*(t-to)+C(t-to)^2
    % O3= A + B*(t-to)+C(t-to)^2+D(t-to)^3
    % to=solar noon
clear all    
%addpath(genpath('c:\Code\iberonesia\Matlab'));
addpath(genpath(fullfile('~','CODE','rbcce.aemet.es','iberonesia','matlab')));
path_root=(fullfile('~','CODE','rbcce.aemet.es','iberonesia','RBCC_E','Triad'))
%warning('off')

%Cargamos los datos/a�os y los agrupamos por brewer en la misma variable.
                  aux157=[];   aux183=[];   aux185=[];
                  aux=cell(3);
for ano=2019:-1:2005
    %�existe el archivo .txt que contiene los datos?
    %a=num2str(ano);
    %a=strcat('c:\Triada\Datos viejos\','triad_comp',a,'.mat')
    %A=exist(a)
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

% Vamos a estudiar dos conjunto de datos para el art�culo: El primero contiene 
% todas las medidas realizadas desde 2005 y en el segundo solo medidas
% sincronizadas. Por tanto, debemos seleccionar de las variables auxXXX
% solo aquellas que se encuentren sincronizadas.

%% 1�Paso: Eliminamos medidas DS con una separaci�n temporal <= 3.293 min.

error157=[];   error183=[];  error185=[];
Tsincro=3.293
% Correcto Tsincro=3.293
for i=1:1:(size(aux157,1)-1)
    if abs(aux157(i,1)-aux157(i+1,1))<=(0.00069444444*Tsincro);
        error157=[error157;i];
    end
end

for i=1:1:(size(aux183,1)-1)
    if abs(aux183(i,1)-aux183(i+1,1))<=(0.00069444444*Tsincro)
        error183=[error183;i];
    end
end

for i=1:1:(size(aux185,1)-1)
    if abs(aux185(i,1)-aux185(i+1,1))<=(0.00069444444*Tsincro)
        error185=[error185;i];
    end
end

error157=[];   error183=[];  error185=[];
% Eliminamos las medidas que no esten muy espaciadas temporalmente.

if isempty(error157)==0 % Por si no hay datos para ese d�a
    aux157(error157(:,1),:)=[];
end
if isempty(error183)==0 % Por si no hay datos para ese d�a
    aux183(error183(:,1),:)=[];
end
if isempty(error185)==0 % Por si no hay datos para ese d�a
    aux185(error185(:,1),:)=[];
end    
% aux183(error183(:,1),:)=[]; aux185(error185(:,1),:)=[];

%Almacenamos las medidas simult�neas. 
aux157s=[];     aux183s=[];       aux185s=[];
Tsincro=3.293
for i=1:1:size(aux183,1)
    t_inicio=aux183(i,1)-0.0017361*Tsincro;
    t_fin=aux183(i,1)+0.0017361*Tsincro;
    a=find(aux185(:,1)>=t_inicio & aux185(:,1)<=t_fin);
    if isempty(a)==0
        b=find(aux157(:,1)>=t_inicio & aux157(:,1)<=t_fin);
        if isempty(b)==0
            aux183s=[aux183s;aux183(i,:)];
            aux157s=[aux157s;aux157(b(1,1),:)];
            aux185s=[aux185s;aux185(a(1,1),:)];
        end        
    end
end

TM=mean([aux157s(:,1),aux183s(:,1),aux185s(:,1)],2);
boxplot(matadd([aux157s(:,1),aux183s(:,1),aux185s(:,1)],-TM)/0.000694444)

clear b; clear error157; clear error183; clear error185; clear i;
clear t_fin; clear t_inicio; clear Tsincro; clear a;

% _______________________________________________________________%

% Hasta aqui tenemos en la variables auxXXX todas las medidas del brewer
% XXX. Mientras que auxXXXs todas las medidas que cumplen la condici�n de
% sincronozaci�n. Es decir, que en 3.293 minutos de separaci�n tengamos las
% medidas de los 3 brewers.

% Las vamos a almacenar en una celda ya que las siguientes operaciones son
% similares para todas. 

Brewer{1}=aux157;  Brewer{2}=aux183;   Brewer{3}=aux185; %Todos las medidas
Brewer{4}=aux157s; Brewer{5}=aux183s;  Brewer{6}=aux185s; % Simult�neas

% Construimos las variables BXXX: A�o (1), d�a Juliano (2), hora de la medida (3)
% ozono (4) y ozono_slcorrec (5), angulo (6),  masa �ptica (7), de las medidas, y una 
% variable de control (8) que permita identificar periodo de tiempos en los 
% que el brewer no estaba Iza�a, averiado, medidas erroneas, etc. 
% La variable de control es 0 por defecto. (0= buen dato, 1= mal dato)

BXXX{1}=[]; BXXX{2}=[]; BXXX{3}=[]; BXXX{4}=[]; BXXX{5}=[]; BXXX{6}=[];
for j=1:1:6
    a=size(Brewer{j})
    BXXX{j}=nan(a(1,1),8);
    %  ---- Brewer 157----  %
    fecha=datevec(Brewer{j}(:,1)); BXXX{j}(:,1)=fecha(:,1); % a�o
    fecha0=fecha; fecha0(:,2)=1; fecha0(:,3)=1;
    fecha0=fix(datenum(fecha0)); BXXX{j}(:,2)=fix(Brewer{j}(:,1)-fecha0+1); % D�a
    BXXX{j}(:,3)=fecha(:,4)*60+fecha(:,5);  % Hora medida (minutos)
    BXXX{j}(:,4)=Brewer{j}(:,4); %ozono 
    BXXX{j}(:,5)=Brewer{j}(:,5); %ozono_slcorrec 
    BXXX{j}(:,6)=Brewer{j}(:,2); %SZA
    BXXX{j}(:,7)=Brewer{j}(:,3); %masa �ptica
    BXXX{j}(:,8)=0; % Variable de control
end
clear fecha; clear fecha0; clear c; clear b; clear a; clear j;

% Definimos las variables BXXX y las sacamos de la celda.
B157=BXXX{1};   B183=BXXX{2};   B185=BXXX{3};
B157s=BXXX{4};  B183s=BXXX{5};  B185s=BXXX{6};

% Almacenamos dentro de la celda "Brewer", los valores de BXXX con el
% objetivo de reducir el tama�o del c�digo de programaci�n ya que el
% proceso de depuraci�n de datos es similar.

Brewer{1}=B157;   Brewer{2}=B183;   Brewer{3}=B185; 
Brewer{4}=B157s;  Brewer{5}=B183s;  Brewer{6}=B185s;

% Depuramos datos: La variable control introducida anteriormente, toma 
% valor distinto de cero si:
  
%   - Si hay menos de 12 medidas diarias y std >0.6 (Hecho) (Control=1)
%   - No hay datos antes y despu�s del medio d�a solar. (Hecho) (Control=2)
%   - Brewer 185 cuando esta de campa�a. (Hecho) (Control=3)
%   - Aver�as claramente identificadas.  (Control=4)

for j=1:1:6
    if j<=3
        inicio=2018; fin=2005;
    else
        inicio=2016; fin=2010;
    end
    for ano=inicio:-1:fin
        Year=find(Brewer{j}(:,1)==ano);
        if isempty(Year)==0  % Para chequear que hay datos para ese a�o
            seleccion=Brewer{j}(Year(1,1):Year(end,1),:);
            for dia=366:-1:1
                day=find(seleccion(:,2)==dia);
                if isempty(day)==0 % Por si no hay datos ese d�a
                    datos=seleccion(day(1,1):day(end,1),:);
                    
                    % --- Filtrado por N�mero de datos y Std ------%
                    if size(datos,1)>=12 && std(datos(:,4))./sqrt(size(datos,1))<=0.6
                        Brewer{j}(Year(1,1)+day(1,1)-1:Year(1,1)+day(end,1)-1,8)=0;
                    else
                        Brewer{j}(Year(1,1)+day(1,1)-1:Year(1,1)+day(end,1)-1,8)=1;
                    end
                    
                    % --- Filtrado: datos antes y despu�s del Solar Noon --- %
                    % Se aplica este SOLO si se pasa el primero.
                    if Brewer{j}(Year(1,1)+day(1,1)-1:Year(1,1)+day(end,1)-1,8)==0;
                        % Solar Noon
                        r_date=datenum(ano,0,0)+dia;
                        noon_m=solar_noon(r_date,-16.499); %Minutos hasta el noon solar
                        auxiliar=seleccion(day(1,1):day(end,1),3)-noon_m;
                        negativos=find(auxiliar<0); % Medidas hechas antes del solar noon.
                        if size(negativos,1)<size(auxiliar,1) && size(negativos,1)>4
                            Brewer{j}(Year(1,1)+day(1,1)-1:Year(1,1)+day(end,1)-1,8)=0;
                        else
                            Brewer{j}(Year(1,1)+day(1,1)-1:Year(1,1)+day(end,1)-1,8)=2;
                        end                      
                    end
                end
            end
        end
    end
end
clear ano; clear auxiliar; clear datos; clear day; clear dia; 
clear negativos; clear noon_m; clear r_date; clear seleccion; clear Year;
clear inicio; clear fin;

% FILTRADO: Brewer 185 en campa�a y aver�as. (Todas y simult�neas)
 
 for j=3:3:6
     % Campa�as
     aa=[];
     a=find(Brewer{j}(:,1)==2007 & Brewer{j}(:,2)>=245 & Brewer{j}(:,2)<=254); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2009 & Brewer{j}(:,2)>=246 & Brewer{j}(:,2)<=265); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2010 & Brewer{j}(:,2)>=200 & Brewer{j}(:,2)<=212); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2011 & Brewer{j}(:,2)>=186 & Brewer{j}(:,2)<=197); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2012 & Brewer{j}(:,2)>=164 & Brewer{j}(:,2)<=181); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2013 & Brewer{j}(:,2)>=155 & Brewer{j}(:,2)<=176); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2014 & Brewer{j}(:,2)>=184 & Brewer{j}(:,2)<=209); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2015 & Brewer{j}(:,2)>=142 & Brewer{j}(:,2)<=157); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2016 & Brewer{j}(:,2)>=175 & Brewer{j}(:,2)<=210); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2017 & Brewer{j}(:,2)>=145 & Brewer{j}(:,2)<=165); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2018 & Brewer{j}(:,2)>=205 & Brewer{j}(:,2)<=225); aa=[aa;a];
     Brewer{j}(aa,8)=3;

     % Aver�as.
     aa=[]
     a=find(Brewer{j}(:,1)==2005 & Brewer{j}(:,2)==276); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2006 & Brewer{j}(:,2)==100); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2006 & Brewer{j}(:,2)==30);  aa=[aa;a];
     a=find(Brewer{j}(:,1)==2006 & Brewer{j}(:,2)==323); aa=[aa;a];
     a=find(Brewer{j}(:,1)==2006 & Brewer{j}(:,2)==75);  aa=[aa;a];
     B185(aa,8)=4;
 end

% FILTRADO: Brewer 183 en campa�a  (Todas y simult�neas)
 
 for j=2:3:5
     % Campa�as
     aa=[];
     a=find(Brewer{j}(:,1)==2013 & Brewer{j}(:,2)>=155 & Brewer{j}(:,2)<=176); aa=[aa;a];
     Brewer{j}(aa,8)=3;
     % Aver�as.
     aa=[];
     a=find(Brewer{j}(:,1)==2006 & Brewer{j}(:,2)>=180& Brewer{j}(:,2)<=282); aa=[aa;a];
     Brewer{j}(aa,8)=4;
 end
 
% Filtrado: Brewer 157 en campa�a y aver�as. (Todas y simult�neas)
for j=1:3:4
    a=find(Brewer{j}(:,1)==2006 & Brewer{j}(:,2)==160); Brewer{j}(a,8)=4;
    a=find(Brewer{j}(:,1)==2014 & Brewer{j}(:,2)>=130 & Brewer{j}(:,2)<=161); 
    Brewer{j}(a,8)=4;
end

% Eliminamos de Brewer{1} todas las medidas descartadas (control >1). 
% Los datos eliminados los almecenamos en la Variable QuitadoXXX. Esto
% �ltimo s�lo para cuando evaluamos TODAS las medidas.

for j=1:1:6
    a=find(Brewer{j}(:,8)>0); 
    if j==1
        Quitado157=Brewer{1}(a,:); save('Quitado157.txt','Quitado157');
    elseif j==2
        Quitado183= Brewer{2}(a,:); save('Quitado183.txt','Quitado183');
    elseif j==3
        Quitado185= Brewer{3}(a,:); save('Quitado185.txt','Quitado185');
    end
    Brewer{j}(a,:)=[];
end

clear Quitado157; clear Quitado183; clear Quitado185; 
clear j; clear a; clear aa;

% Reescribimos las variables iniciales BXXX (Todas y simult�neas).

B157=Brewer{1};   B183=Brewer{2};   B185=Brewer{3}; 
B157s=Brewer{4};  B183s=Brewer{5};  B185s=Brewer{6}; 

%% Comparativa m�todo triada C�nada (Polinomio de 2� y 3� grado)

% Daily model for midday
% O3= A + B*(t-to) + C(t-to)^2
% to=solar noon

% Variables donde almacenamos los datos mediaXXX:=[a�o,dia,A(2�,O3),
% Asl(2�,O3sl),A(3�,O3),Asl(3�,O3sl),media(O3),media(O3sl),control]
% Introducimos una variable control que indique cuando debemos emplear 
% los coeficientes calculados a partir del ozono corregido por SL. En este
% apartado la variable control toma el valor cero siempre (control=0). Ser�
% en la siguiente secci�n donde se a�adida los casos donde Control =1. Es
% decir, los periodos donde deberems coger los coeficientes calculados
% teniendo en cuenta la correcci�n por Sl.

mediaXXX{1}=[]; mediaXXX{2}=[]; mediaXXX{3}=[]; 
mediaXXX{4}=[]; mediaXXX{5}=[]; mediaXXX{6}=[];  

for j=1:1:6
    if j<=3
        inicio=2018; fin=2005;
    else
        inicio=2016; fin=2010;
    end
    for ano=inicio:-1:fin
        for dia=366:-1:1
            a=find(Brewer{j}(:,1)==ano & Brewer{j}(:,2)==dia); % selecci�n de a�o y d�a
            if isempty(a)==0 % Por si no hay datos para ese d�a
                datos=Brewer{j}(a,:); % Datos que cumplen la condici�n del "find"
                % Solar Noon.
                r_date=datenum(ano,0,0)+dia;
                noon=solar_noon(r_date, -16.499);
                idx=ones(size(a,1),1);
                % Medidas ozono
                d2=[datos(:,4),idx,datos(:,3)-noon,(datos(:,3)-noon).^2]; %[ozono,(t-t0),(t-t0)^2] 
                d3=[datos(:,4),idx,datos(:,3)-noon,(datos(:,3)-noon).^2,(datos(:,3)-noon).^3]; %[ozono,(t-t0),(t-t0)^2,(t-t0)^3]
                [B2]=regress(d2(:,1),d2(:,2:end));  % Ajuste polinomio grado 2
                [B3]=regress(d3(:,1),d3(:,2:end));  % Ajuste polinomio grado 3
                % Medidas ozono_slcorrec
                d2sl=[datos(:,5),idx,datos(:,3)-noon,(datos(:,3)-noon).^2]; %[ozono,(t-t0),(t-t0)^2] 
                d3sl=[datos(:,5),idx,datos(:,3)-noon,(datos(:,3)-noon).^2,(datos(:,3)-noon).^3]; %[ozono,(t-t0),(t-t0)^2,(t-t0)^3]
                [B2sl]=regress(d2sl(:,1),d2(:,2:end));  % Ajuste polinomio grado 2
                [B3sl]=regress(d3sl(:,1),d3(:,2:end));  % Ajuste polinomio grado 3
                mediaXXX{j}=[mediaXXX{j};ano,dia,B2(1,1),B2sl(1,1),B3(1,1),B3sl(1,1),mean(datos(:,4)),mean(datos(:,5)),0];           
        else % Si no hay datos para este d�a
            aux1=nan(1,6);
            mediaXXX{j}=[mediaXXX{j};ano,dia,aux1,0]; 
        end
    end
end
end

%% Modificamos la variable control A

% Nota: En las variables mediaXXX se encuentra los valores del coefficiente A 
% del ajuste de fioletov (polinomio 2 y 3 grado). Asi como el valor medio.
% No obstante, estas cuentas han sido hechas para los valores de ozono y 
% ozono_slcorrec. Por tanto, podr�amos hablar de coeficientes A(ozono) y 
% Asl(ozono_slcorrec).

% En la ultima columna de mediaXXX esta la variable de control, que por 
% defecto toma valor cero (control=0). Esto indica que se seleccionar�
% siempre los valores calculados sin la correcci�n por SL. 

% En esta secci�n, vamos a a�adir los periodos de tiempo en que deben de
% a�adirse las correcciones por SL. (variable control =1) frente a los
% valores sin correcci�n por SL.

for j=1:1:6
   if j==1 
   a=find(mediaXXX{j}(:,1)==2014 & mediaXXX{j}(:,2)>=258 & mediaXXX{j}(:,2)<=366);
   mediaXXX{j}(a,9)=1;
   a=find( mediaXXX{j}(:,1)==2016 &  mediaXXX{j}(:,2)>=1 &  mediaXXX{j}(:,2)<=99);
   mediaXXX{j}(a,9)=1;
   elseif j==4
       a=find(mediaXXX{j}(:,1)==2014 & mediaXXX{j}(:,2)>=258 & mediaXXX{j}(:,2)<=366);
       mediaXXX{j}(a,9)=1;
       a=find( mediaXXX{j}(:,1)==2016 &  mediaXXX{j}(:,2)>=1 &  mediaXXX{j}(:,2)<=99);
       mediaXXX{j}(a,9)=1;
   elseif j==2 
       a=find(mediaXXX{j}(:,1)==2006 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2007 & mediaXXX{j}(:,2)>=179 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2008 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2009 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=232);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2012 & mediaXXX{j}(:,2)>=345 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2013 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2014 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=335); mediaXXX{j}(a,9)=1;
    elseif j==5 
       aa=[];
       a=find(mediaXXX{j}(:,1)==2006 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2007 & mediaXXX{j}(:,2)>=179 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2008 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=366); mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2009 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=232);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2012 & mediaXXX{j}(:,2)>=345 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2013 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2014 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=335); mediaXXX{j}(a,9)=1;    
   elseif j==3 
       a=find(mediaXXX{j}(:,1)==2011 & mediaXXX{j}(:,2)>=312 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2012 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=61);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2014 & mediaXXX{j}(:,2)>=44 & mediaXXX{j}(:,2)<=343);  mediaXXX{j}(a,9)=1;
   elseif j==6
       a=find(mediaXXX{j}(:,1)==2011 & mediaXXX{j}(:,2)>=312 & mediaXXX{j}(:,2)<=366);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2012 & mediaXXX{j}(:,2)>=1 & mediaXXX{j}(:,2)<=61);  mediaXXX{j}(a,9)=1;
       a=find(mediaXXX{j}(:,1)==2014 & mediaXXX{j}(:,2)>=44 & mediaXXX{j}(:,2)<=343);  mediaXXX{j}(a,9)=1;
   end
end
clear a; clear aa; clear ano; clear ans; clear aux1; clear B2; clear B2sl;
clear B3; clear B3sl; clear d2; clear d2sl; clear d3; clear d3sl; clear datos;
clear dia; clear fin; clear inicio; clear j; clear idx; clear noon;
clear r_date;

% Selecionamos los coeficientes A en funci�n de la variable control. Es
% decir, si cogemos los datos con SL o sin SL correction.
% media{}=[a�o, d�a, A2, A3, media].

aux=mediaXXX;
media{1}=[]; media{2}=[]; media{3}=[]; media{4}=[]; media{5}=[]; media{6}=[];

for j=1:1:6
    media{j}=nan(size(aux{j},1),5);
    a=find(aux{j}(:,9)==0);
    media{j}(a,:)=[aux{j}(a,1:2),aux{j}(a,3),aux{j}(a,5),aux{j}(a,7)]; % Datos sin SL correct
    a=find(mediaXXX{j}(:,9)==1);
    media{j}(a,:)=[aux{j}(a,1:2),aux{j}(a,4),aux{j}(a,6),aux{j}(a,8)]; % Datos con SL correct
end

clear j; clear aux; clear a;
save workspace_5min_report
%% Pintamos los resultados del ajuste de Fieletov

% HISTROGRAMA ERROR DIARIO (solo empleamos los datos cuando los tres brewers han medido ese dia)

[Histograma,error,error_relativo]=histog(media,1500,600,550,550)

% HISTROGRAMA ERROR MENSUAL

%Todas las medidas
mes_A2=mensual_A2(media)
mes_A3=mensual_A3(media)
mes_AM=mensual_AM(media)

%Simultaneas
mes_A2s=mensual_A2s(media)
mes_A3s=mensual_A3s(media)
mes_AMs=mensual_AMs(media)

[error,desviacion]=histog_mensual(mes_A2,mes_A3,mes_AM)
[error_simultaneas,std_simultaneas]=histog_mensual(mes_A2s,mes_A3s,mes_AMs)


% Figura 2.

% Cogemos los datos de mes_AM (promedio mensual de las medias diarias)

figure
ejex=datenum(media{1}(:,1),1,media{1}(:,2))
plot(ejex,media{1}(:,3),'o'); hold on

aux=[]
for j=1:1:size(mes_AM,1)
    if isnan(mes_AM(j,3))==0
        aux=[aux;mes_AM(j,:)]
    end
end
 
ejex=datenum(aux(:,1),aux(:,2),1)
plot(ejex,aux(:,3),'LineWidth',3)

%legend('BR#157','BR#183','BR#185');
ylabel('Total ozone (DU)')
xlabel('Year')
ylim([230 350])
startDate = datenum('01-07-2005');
endDate = datenum('01-01-2019');
xData = linspace(startDate,endDate,84);
datetick('x','yyyy','keeplimits')
xlim([datenum(2005,7,1),datenum(2019,1,1)])

a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',12);


%% Grafica report Iza�a

% Funci�n empleada para el ajuste Histograma y funci�n gaussiana. 
% Futuro Paper de la Triada.

%Par�metros entrada
    
    % media=celda que contiene, media{}=[a�o, d�a, A2, A3, media].
    % para todas las medidas y las simult�neas.
    
    % EjeX1 y EjeX2 son los valores de corte para los histogramas

periodo=[1, 31; 32, 59; 60, 90; 91, 120; 121,151; 152,181; 182,212; 213,243; 244, 273;274, 304; 305, 334; 335, 366];

medias=[media{1}(:,1:3),media{2}(:,3),media{3}(:,3)]

mensual=[]
contador=[];
for ano=2019:-1:2005
    for j=12:-1:1
         a=find(medias(:,1)==ano & medias(:,2)>=periodo(j,1) & medias(:,2)<=periodo(j,2))
         datos=[medias(a,:),sum(isnan(medias(a,:)),2)]
         a=find(datos(:,6)==0)
         contador=[size(a,1),contador]
         if size(a,1)>=3
            seleccinados=datos(a,:) %D�as y ozono en que los brewers tiene medidas simult�neas.
            % Calculamos el valor medio Me
            Me157=[ano,j,mean(seleccinados(:,3))]; 
            Me183=[ano,j,mean(seleccinados(:,4))];
            Me185=[ano,j,mean(seleccinados(:,5))]; 
            b=[ano,j,Me157(1,3),Me183(1,3),Me185(1,3)]; mensual=[mensual;b]
         else
            b=[ano,j,nan,nan,nan]; mensual=[mensual;b] 
         end 
    end
end

% Tenemos en cuenta los meses donde el brewer 183 estuvo roto (Delta)

a=find(medias(:,1)==2005 & medias(:,2)>=periodo(12,1) & medias(:,2)<=periodo(12,2))
b1=[medias(a,1:3),medias(a,5)]
datos=[b1,sum(isnan(b1),2)]
a=find(datos(:,5)==0)
if size(a,1)>=3
    seleccinados=datos(a,:) %D�as y ozono en que los brewers tiene medidas simult�neas.
    Me157=mean(seleccinados(:,3));
    Me185=mean(seleccinados(:,4));
    a=find(mensual(:,1)==2005 & mensual(:,2)==12);
    mensual(a,3)=Me157;  mensual(a,5)=Me185;
end

for j=1:1:7
    a=find(medias(:,1)==2006 & medias(:,2)>=periodo(j,1) & medias(:,2)<=periodo(j,2))
    b1=[medias(a,1:3),medias(a,5)]
    datos=[b1,sum(isnan(b1),2)]
    a=find(datos(:,5)==0)
    if size(a,1)>=3
        seleccinados=datos(a,:) %D�as y ozono en que los brewers tiene medidas simult�neas.
        Me157=mean(seleccinados(:,3));
        Me185=mean(seleccinados(:,4));
        a=find(mensual(:,1)==2006 & mensual(:,2)==j);
        mensual(a,3)=Me157;  mensual(a,5)=Me185;
    end
end

Emensual=[]
for ano=2018:-1:2005
    for j=12:-1:1
        a=find(mensual(:,1)==ano & mensual(:,2)==j)
        if isempty(a)==0
            c=[mensual(a,3),mensual(a,4),mensual(a,5)]
            media1=nanmean(c)
            Emensual=[ano,j,mensual(a,3),mensual(a,4),mensual(a,5),media1,(mensual(a,3)/media1),(mensual(a,4)/media1),(mensual(a,5)/media1);Emensual];
        end
    end
end

%mensual(85:144,:)=[]
%Emensual(1:60,:)=[]
%save Emensual.mat Emensual
figure
subplot(2,1,1)
ejex=datenum(mensual(:,1),1,1);
ejex=datenum(mensual(:,1),mensual(:,2),1);
plot(ejex,mensual(:,3),'o-','linewidth',2); hold on
plot(ejex,mensual(:,4),'*-','linewidth',2); hold on
plot(ejex,mensual(:,5),'*-','linewidth',2); hold on
legend('BR#157','BR#183','BR#185');
ylabel('Total ozone (DU)')
ylim([240 340])
startDate = datenum('01-07-2005');
endDate = datenum('01-01-2019');
xData = linspace(startDate,endDate,84);
datetick('x','yyyy','keeplimits')
xlim([datenum(2005,7,1),datenum(2019,1,1)])


%title(sprintf('O3 Monthly in the Period 2005-2018'));
subplot(2,1,2)
ejex(:,1)=datenum(Emensual(:,1),Emensual(:,2),1);
plot(ejex,Emensual(:,7),'+-',ejex,Emensual(:,8),'o-',ejex,Emensual(:,9),'*-','linewidth',2);
legend('BR#157','BR#183','BR#185');
startDate = datenum('01-07-2005');
endDate = datenum('01-01-2019');
xData = linspace(startDate,endDate,84);
datetick('x','yyyy','keeplimits')
ylabel('Deviation (%)')
% x=ejex;
% datetick('x','yyyy','Keepticks','keeplimits');
% xlim([datenum(2010,1,1),datenum(2016,12,31)])
% legend('BR#157','BR#183','BR#185');
xlim([datenum(2005,7,1),datenum(2019,1,1)])
%title(sprintf('O3 Error Monthly in the Period 2005-2015'));

a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',12);


figure
datos_Histograma=Histograma{2}(:,6:8)
nbins = 20;
hist(datos_Histograma,nbins)
%title(sprintf('A3  O3 Daily Error respect to mean value in the Period 2005-2018'))
legend('BR#157','BR#183','BR#185');
ylabel('Number of observations (Days)')
xlabel('Daily difference (DU)')
xlim([-5 5])

a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',12);

% 
% 
% a=findobj(gcf);
% 
% 
% figure
% subplot(1,2,1)
% datos_Histograma=Histograma{1}(:,6:8)
% nbins = 14;
% [count, centers]=hist(datos_Histograma,nbins)
% %title(sprintf('Daily Error respect to mean value in the Period 2005-2016'))
% title(sprintf('Difference respect to mean value in the Period 2010-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations (Days)')
% xlabel('Average daily difference (DU')
% xlim([-4 4])
% 
%  % get the handles associated with the current figure
% 
% 
% subplot(1,2,2)
% datos_Histograma=Histograma{2}(:,6:8)
% nbins = 20;
% hist(datos_Histograma,nbins)
% title(sprintf('A3  O3 Daily Error respect to mean value in the Period 2005-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations (Days)')
% xlim([-5 5])
% % Todas las medidas y las simult�neas!!! Coeficiente A2 y A3
% % Periodo 2010-2016
% clear aux
% clear Histograma
% aux{1}=[media{1}(1:2196,1:3),media{2}(1:2196,3),media{3}(1:2196,3)]% A2, Todas medidas.
% aux{2}=[media{1}(1:2196,1:2),media{1}(1:2196,4),media{2}(1:2196,4),media{3}(1:2196,4)]% A3, Todas medidas.
% aux{3}=[media{1}(1:2196,1:2),media{1}(1:2196,5),media{2}(1:2196,5),media{3}(1:2196,5)]% Media, Todas medidas.
% 
% aux{4}=[media{4}(1:2196,1:3),media{5}(1:2196,3),media{6}(1:2196,3)]% A2, Simult�neas.
% aux{5}=[media{4}(1:2196,1:2),media{4}(1:2196,4),media{5}(1:2196,4),media{6}(1:2196,4)]% A3, Simult�neas.
% aux{6}=[media{4}(1:2196,1:2),media{4}(1:2196,5),media{5}(1:2196,5),media{6}(1:2196,5)]% Media, Simult�neas.
% 
% Histograma{1}=[]; Histograma{2}=[]; Histograma{3}=[];
% Histograma{4}=[]; Histograma{5}=[]; Histograma{6}=[];
% 
% for j=1:1:6
%     b=[]
%     for i=1:1:size(aux{j},1)
%         if sum(isnan(aux{j}(i,:)))==0
%             a=mean(aux{j}(i,3:5))
%             if abs(aux{j}(i,3)-a)<=5 && abs(aux{j}(i,4)-a)<=5 && abs(aux{j}(i,5)-a)<=5
%                 b=[aux{j}(i,1:2),aux{j}(i,3),aux{j}(i,4),aux{j}(i,5),aux{j}(i,3)-a,aux{j}(i,4)-a,aux{j}(i,5)-a;b];
%             end
%         end
%     end
%     Histograma{j}=b   
% end
% 
% figure
% subplot(3,3,[1 4])
% datos_Histograma=Histograma{1}(:,6:8)
% nbins = 20;
% hist(datos_Histograma,nbins)
% title(sprintf('A2 Todas medidas  Period 2010-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations')
% xlim([-5 5])
% 
% subplot(3,3,7)
% datos_Histograma=Histograma{2}(:,6:8)
% nbins = 20;
% hist(datos_Histograma,nbins)
% title(sprintf('A3  Todas medidas Period 2010-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations')
% xlim([-5 5])
% 
% subplot(2,3,4)
% datos_Histograma=Histograma{4}(:,6:8)
% nbins = 20;
% hist(datos_Histograma,nbins)
% title(sprintf('A2 Simult�neas Period 2010-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations')
% xlim([-5 5])
% 
% subplot(2,3,5)
% datos_Histograma=Histograma{5}(:,6:8)
% nbins = 20;
% hist(datos_Histograma,nbins)
% title(sprintf('A3 Simult�neas  Period 2010-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations')
% xlim([-5 5])
% 
% subplot(2,3,3)
% datos_Histograma=Histograma{3}(:,6:8)
% nbins = 20;
% hist(datos_Histograma,nbins)
% title(sprintf('Media Period 2010-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations')
% xlim([-5 5])
% 
% subplot(2,3,6)
% datos_Histograma=Histograma{6}(:,6:8)
% nbins = 20;
% hist(datos_Histograma,nbins)
% title(sprintf('Media Simult�neas Period 2010-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations')
% xlim([-5 5])
% 
% 
% %% Grafica mensual
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
% 
% %% Busqueda de errores sistem�ticos. 
% 
% aux=[media{1}(:,1),media{1}(:,2),media{1}(:,5),media{2}(:,5),media{3}(:,5)] % Cojo los valores medios.
% datos=[];
% 
% for i=1:1:size(aux,1)
%     if sum(isnan(aux(i,:)))==0
%         datos=[datos;aux(i,:)];
%     end
% end
% 
% orden=ones(size(datos,1),3);
% for i=1:1:size(datos,1)
%     aa=sort(datos(i,3:5));
%     a1=find(aa==datos(i,3)); orden(i,a1)=157;
%     a2=find(aa==datos(i,4)); orden(i,a2)=183;
%     a3=find(aa==datos(i,5)); orden(i,a3)=185;  
% end
% 
% Repeticion=ones(1,150)
% 
% referencia=orden(1,:)
% cont=1
% for i=2:1:size(datos,1)
%     if referencia==orden(i,:)
%         cont=cont+1
%     end
%     if referencia~=orden(i,:)
%         Repeticion(1,cont)=Repeticion(1,cont)+1
%         referencia=orden(i,:)
%         cont=1
%     end
% end
% 
% Repeticion=Repeticion-ones(1,150)
% total=0
% for i=1:1:150
%     total1=i*Repeticion(1,i)
%     total=total+total1
% end
% 
% plot([1:1:150], Repeticion')
% nbins = 200;
% hist(Repeticion',nbins)
% sum(Repeticion)
% xlim([0 10])
% aa=sort(datos(1,3:5))
% %% Datos para Omaira y comparativa entre ALberto, eubrewnet y mis datos
% 
% % En esta secci�n preparo los datos que le voy a dar a Omaira. Para ello,
% % al comienzo del codigo guarde las variables auxXXXoma. 
% 
% B{1}=aux157oma; B{2}=aux183oma; B{3}=aux185oma;
% BXXXoma{1}=[]; BXXXoma{2}=[]; BXXXoma{3}=[]; 
% 
% for j=1:1:3
%     a=size(B{j});
%     BXXXoma{j}=nan(a(1,1),8);
%     
%     fecha=datevec(B{j}(:,1)); BXXXoma{j}(:,1)=fecha(:,1); % a�o
%     fecha0=fecha; fecha0(:,2)=1; fecha0(:,3)=1;
%     fecha0=fix(datenum(fecha0)); BXXXoma{j}(:,2)=fix(B{j}(:,1)-fecha0+1); % D�a
%     BXXXoma{j}(:,3)=fecha(:,4)*60+fecha(:,5);  % Hora medida (minutos)
%     BXXXoma{j}(:,4)=B{j}(:,4); %ozono 
%     BXXXoma{j}(:,5)=B{j}(:,5); %ozono_slcorrec 
%     BXXXoma{j}(:,6)=B{j}(:,2); %SZA
%     BXXXoma{j}(:,7)=B{j}(:,3); %masa �ptica
%     BXXXoma{j}(:,8)=0; % Variable de control
% end
% 
% aux157oma=BXXXoma{1}; aux183oma=BXXXoma{2}; aux185oma=BXXXoma{3};
% % FILTRADO: Eliminamos las medidas realizadas en campa�as de calibraci�n 
% % o cuanto el equipo esta averiado.
% 
% %Brewer#185
% 
% aa=[];
% a=find(aux185oma(:,1)==2007 & aux185oma(:,2)>=245 & aux185oma(:,2)<=254); aa=[aa;a];
% a=find(aux185oma(:,1)==2009 & aux185oma(:,2)>=246 & aux185oma(:,2)<=265); aa=[aa;a];
% a=find(aux185oma(:,1)==2010 & aux185oma(:,2)>=200 & aux185oma(:,2)<=212); aa=[aa;a];
% a=find(aux185oma(:,1)==2011 & aux185oma(:,2)>=186 & aux185oma(:,2)<=197); aa=[aa;a];
% a=find(aux185oma(:,1)==2012 & aux185oma(:,2)>=164 & aux185oma(:,2)<=181); aa=[aa;a];
% a=find(aux185oma(:,1)==2013 & aux185oma(:,2)>=155 & aux185oma(:,2)<=176); aa=[aa;a];
% a=find(aux185oma(:,1)==2014 & aux185oma(:,2)>=184 & aux185oma(:,2)<=209); aa=[aa;a];
% a=find(aux185oma(:,1)==2016 & aux185oma(:,2)>=142 & aux185oma(:,2)<=157); aa=[aa;a];      
% a=find(aux185oma(:,1)==2005 & aux185oma(:,2)==276); aa=[aa;a];
% a=find(aux185oma(:,1)==2006 & aux185oma(:,2)==100); aa=[aa;a];
% a=find(aux185oma(:,1)==2006 & aux185oma(:,2)==30);  aa=[aa;a];
% a=find(aux185oma(:,1)==2006 & aux185oma(:,2)==323); aa=[aa;a];
% a=find(aux185oma(:,1)==2006 & aux185oma(:,2)==75);  aa=[aa;a];
% aux185oma(aa,:)=[];
% B{1}(aa,:)=[];
% %Brewer#183
% aa=[];
% a=find(aux183oma(:,1)==2013 & aux183oma(:,2)>=155 & aux183oma(:,2)<=176); aa=[aa;a];       
% a=find(aux183oma(:,1)==2006 & aux183oma(:,2)>=180& aux183oma(:,2)<=282); aa=[aa;a]; 
% aux183oma(aa,:)=[];
% B{2}(aa,:)=[]
% %Brewer#157
% aa=[];
% a=find(aux157oma(:,1)==2006 & aux157oma(:,2)==160); 
% aux157oma(a,:)=[]; B{3}(a,:)=[]
% a=find(aux157oma(:,1)==2014 & aux157oma(:,2)>=130 & aux157oma(:,2)<=161); 
% aux157oma(a,:)=[]; B{3}(a,:)=[]
% 
% % Debemos seleccionar los periodos donde se emplea el ozono o ozono_sl correc
% B{1}=[B{1},(ones(size(B{1},1),1)-1)]
% B{2}=[B{2},(ones(size(B{2},1),1)-1)]
% B{3}=[B{3},(ones(size(B{3},1),1)-1)]
% 
% % Brewer 157
% a=find(aux157oma(:,1)==2014 & aux157oma(:,2)>=258 & aux157oma(:,2)<=366);
% aux157oma(a,8)=1; B{1}(a,6)=1;
% a=find(aux157oma(:,1)==2016 & aux157oma(:,2)>=1 & aux157oma(:,2)<=99);
% aux157oma(a,8)=1; B{1}(a,6)=1;
% 
% % Brewer 183 
% a=find(aux183oma(:,1)==2006 & aux183oma(:,2)>=1 & aux183oma(:,2)<=366);  aux183oma(a,8)=1; B{2}(a,6)=1;
% a=find(aux183oma(:,1)==2007 & aux183oma(:,2)>=179 & aux183oma(:,2)<=366); aux183oma(a,8)=1; B{2}(a,6)=1; 
% a=find(aux183oma(:,1)==2008 & aux183oma(:,2)>=1 & aux183oma(:,2)<=366); aux183oma(a,8)=1; B{2}(a,6)=1;
% a=find(aux183oma(:,1)==2009 & aux183oma(:,2)>=1 & aux183oma(:,2)<=232); aux183oma(a,8)=1; B{2}(a,6)=1;
% a=find(aux183oma(:,1)==2012 & aux183oma(:,2)>=345 & aux183oma(:,2)<=366); aux183oma(a,8)=1; B{2}(a,6)=1;
% a=find(aux183oma(:,1)==2013 & aux183oma(:,2)>=1 & aux183oma(:,2)<=366); aux183oma(a,8)=1; B{2}(a,6)=1;
% a=find(aux183oma(:,1)==2014 & aux183oma(:,2)>=1 & aux183oma(:,2)<=335); aux183oma(a,8)=1; B{2}(a,6)=1;
%     
% % Brewer 185
% a=find(aux185oma(:,1)==2011 & aux185oma(:,2)>=312 & aux185oma(:,2)<=366); aux185oma(a,8)=1; B{3}(a,6)=1;
% a=find(aux185oma(:,1)==2012 & aux185oma(:,2)>=1 & aux185oma(:,2)<=61); aux185oma(a,8)=1; B{3}(a,6)=1;
% a=find(aux185oma(:,1)==2014 & aux185oma(:,2)>=44 & aux185oma(:,2)<=343); aux185oma(a,8)=1; B{3}(a,6)=1;
% 
% %  Guardamos los datos para Omaira. 
% save('ozono157sergio.mat','aux157oma');
% save('ozono183sergio.mat','aux183oma');
% save('ozono185sergio.mat','aux185oma');
% 
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