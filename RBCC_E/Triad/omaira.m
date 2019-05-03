function [error,desviacion]=Rene(Brewer)


% En la matriz Brewer, existen dos columnas de ozono, uno calculado con la
% SL correction y otros que no. Debemos seleccionar el ozono correcto en
% cada caso, por equipo.

% Creamos una matriz Brewer_corregido. Donde solo aparece una única columna
% de ozono, que será la que emplearemos para hacer los calculos.

for j=1:3:4 %Brewer 157, todas las medidas y simultanes
    a=find(Brewer{j}(:,1)==2014 & Brewer{j}(:,2)>=258 & Brewer{j}(:,2)<=366);
    Brewer{j}(a,9)=1;
    a=find( Brewer{j}(:,1)==2016 &  Brewer{j}(:,2)>=1 &  Brewer{j}(:,2)<=99);
    Brewer{j}(a,9)=1;
end

for j=2:3:5 % Brewer 183, todas las medidas y simultanes
    a=find(Brewer{j}(:,1)==2006 & Brewer{j}(:,2)>=1 & Brewer{j}(:,2)<=366);  Brewer{j}(a,9)=1;
    a=find(Brewer{j}(:,1)==2007 & Brewer{j}(:,2)>=179 & Brewer{j}(:,2)<=366);  Brewer{j}(a,9)=1;
    a=find(Brewer{j}(:,1)==2008 & Brewer{j}(:,2)>=1 & Brewer{j}(:,2)<=366);  Brewer{j}(a,9)=1;
    a=find(Brewer{j}(:,1)==2009 & Brewer{j}(:,2)>=1 & Brewer{j}(:,2)<=232);  Brewer{j}(a,9)=1;
    a=find(Brewer{j}(:,1)==2012 & Brewer{j}(:,2)>=345 & Brewer{j}(:,2)<=366);  Brewer{j}(a,9)=1;
    a=find(Brewer{j}(:,1)==2013 & Brewer{j}(:,2)>=1 & Brewer{j}(:,2)<=366);  Brewer{j}(a,9)=1;
    a=find(Brewer{j}(:,1)==2014 & Brewer{j}(:,2)>=1 & Brewer{j}(:,2)<=335); Brewer{j}(a,9)=1;
end

for j=3:3:6 % Brewer 185, todas las medidas y simultanes
    a=find(Brewer{j}(:,1)==2011 & Brewer{j}(:,2)>=312 & Brewer{j}(:,2)<=366);  Brewer{j}(a,9)=1;
    a=find(Brewer{j}(:,1)==2012 & Brewer{j}(:,2)>=1 & Brewer{j}(:,2)<=61);  Brewer{j}(a,9)=1;
    a=find(Brewer{j}(:,1)==2014 & Brewer{j}(:,2)>=44 & Brewer{j}(:,2)<=343);  Brewer{j}(a,9)=1;
end

for j=1:1:6
    Brewer_corr{j}=nan(size(Brewer{j},1),4);
    a=find(Brewer{j}(:,9)==0);
    Brewer_corr{j}(a,:)=[Brewer{j}(a,1:4)]; % Datos sin SL correct
    a=find(Brewer{j}(:,9)==1);
    Brewer_corr{j}(a,:)=[Brewer{j}(a,1:3),Brewer{j}(a,5)]; % Datos con SL correct
end

datos=[Brewer_corr{1};Brewer_corr{2};Brewer_corr{3}]
ajuste2=[];  ajuste3=[]
for k=2016:-1:2005
    for i=366:-1:1
       a=find(datos(:,1)==k & datos(:,2)==i)
       if isempty(a)==0
       % Solar Noon.
       r_date=datenum(k,0,0)+i;
       noon=solar_noon(r_date, -16.499);
       idx=ones(size(a,1),1);
       % Medidas ozono
       d2=[datos(a,4),idx,datos(a,3)-noon,(datos(a,3)-noon).^2]; %[ozono,(t-t0),(t-t0)^2] 
       d3=[datos(a,4),idx,datos(a,3)-noon,(datos(a,3)-noon).^2,(datos(a,3)-noon).^3]; %[ozono,(t-t0),(t-t0)^2,(t-t0)^3]
       [B2]=regress(d2(:,1),d2(:,2:end));  % Ajuste polinomio grado 2
       [B3]=regress(d3(:,1),d3(:,2:end));  % Ajuste polinomio grado 3
       ajuste2=[k,i,B2';ajuste2]; ajuste3=[k,i,B3';ajuste3] 
       end
    end
end


% Calculo el desplazamiento de cada brewer con respecto al valor de la
% triada
delta=[]
desv=[]
delta_relativo=[]
desv_relativa=[]
for i=1:1:size(ajuste2,1)
    ano=ajuste2(i,1);
    dia=ajuste2(i,2);
    for j=1:1:3
        a=find(Brewer_corr{j}(:,1)==ano & Brewer_corr{j}(:,2)==i)
        if isempty(a)==0 % Si tenemos datos
            r_date=datenum(ano,0,0)+dia; 
            noon=solar_noon(r_date, -16.499)
            time=(Brewer{j}(a,3)-noon)
            X=fliplr(ajuste2(i,3:end))
            O3_teorico=polyval(X,time)
            valores=(Brewer{j}(a,4)-O3_teorico)
            mean(valores)
            O3_teorico2(:,1)=ajuste2(:,3)+ajuste2(:,3)*
        
    
    
    
    
for k=2016:-1:2005
    for i=366:-1:1
       Valores2=[]; Valores3=[]
       a=find(ajust2(:,1)==k & ajuste2(:,2)==i)
       if isempty(a)==0
       r_date=datenum(k,0,0)+i;
       noon=solar_noon(r_date, -16.499)
       Valores2(:,1)=ajuste2(:,3)+ajuste2(:,3)*
       
       
       













        