function [mes_A2s]=mensual_A2s(media)

% Funci�n empleada para el ajuste Histograma y funci�n gaussiana. 
% Futuro Paper de la Triada.

%Par�metros entrada
    
    % media=celda que contiene, media{}=[a�o, d�a, A2, A3, media].
    % para todas las medidas y las simult�neas.
    
    % EjeX1 y EjeX2 son los valores de corte para los histogramas

periodo=[1, 31; 32, 59; 60, 90; 91, 120; 121,151; 152,181; 182,212; 213,243; 244, 273;274, 304; 305, 334; 335, 366];

medias=[media{4}(:,1:2),media{4}(:,4),media{5}(:,4),media{6}(:,4)]

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
for ano=2019:-1:2005
    for j=12:-1:1
        a=find(mensual(:,1)==ano & mensual(:,2)==j)
        if isempty(a)==0
            c=[mensual(a,3),mensual(a,4),mensual(a,5)]
            media1=nanmean(c)
            Emensual=[ano,j,mensual(a,3),mensual(a,4),mensual(a,5),media1,(mensual(a,3)/media1),(mensual(a,4)/media1),(mensual(a,5)/media1);Emensual];
        end
    end
end
mes_A2s=Emensual
end