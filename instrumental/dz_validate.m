clear all
run(fullfile('..','read_config_'))

brewer=[157,183,185];
brwid=3

%Date.dayend=now; Date.day0=datenum(2019,1,1); % Esto para los weekly
Date.day0=datenum(2019,2,1);
Date.dayend=datenum(2019,3,1); 

Date.CALC_DAYS =Date.day0:Date.dayend;
Cal.Date=Date;

% READ Bfiles dz
ozone{brwid}=read_bdata_dz(brwid,Cal);

t=[]; a=[]; b=[]; ab=[]; cy=[]; z=[];
nd=length(ozone{brwid}.dz_raw0)

%leemos la slit 2 en a, slit 5 en b y la slit 7 (doble) en ab, dark en z, el número de ciclos en cy, y la fecha en t
for i=1:nd
    t=[t;ozone{brwid}.dz_raw0{i,1}(:,1)];
    cy=[cy;ozone{brwid}.dz_raw0{i,1}(:,10)];
    z=[z;ozone{brwid}.dz_raw0{i,1}(:,12)];
    a=[a;ozone{brwid}.dz_raw0{i,1}(:,14)];
    b=[b;ozone{brwid}.dz_raw0{i,1}(:,16)];
    ab=[ab;ozone{brwid}.dz_raw0{i,1}(:,18)];
end

% pasamos a count rates
it=0.1147;
data=[t,cy,z,a,b,ab];
data=[data 2*(a-z)./cy/it 2*(b-z)./cy/it 2*(ab-z)./cy/it];

% leemos el dt para cada fecha de las configuraciones y lo ponemos en las
% columas 10 (op) y 11 (alt)
for i=1:length(t)
    if t(i)>icf_op{brwid}(2,end)
        data(i,10)=icf_op{brwid}(14,end);
    else
        data(i,10)=icf_op{brwid}(14,find(icf_op{brwid}(2,:)>t(i),1)-1);
    end
    if t(i)>icf_a{brwid}(2,end)
        data(i,10)=icf_a{brwid}(14,end);
    else
        data(i,11)=icf_a{brwid}(14,find(icf_a{brwid}(2,:)>t(i),1)-1);
    end
end

%corregimos las medidas con el deadtine operativo
data(:,12)=data(:,7);
data(:,13)=data(:,8);
data(:,14)=data(:,9);
for i=1:9
    data(:,12)=data(:,7).*exp(data(:,12)*data(i,10));
    data(:,13)=data(:,8).*exp(data(:,13)*data(i,10));
    data(:,14)=data(:,9).*exp(data(:,14)*data(i,10));
end
data(:,15)=(data(:,12)+data(:,13))./data(:,14);

%corregimos las medidas con el deadtine alternativo
data(:,16)=data(:,7);
data(:,17)=data(:,8);
data(:,18)=data(:,9);
for i=1:9
    data(:,16)=data(:,7).*exp(data(:,16)*data(i,11));
    data(:,17)=data(:,8).*exp(data(:,17)*data(i,11));
    data(:,18)=data(:,9).*exp(data(:,18)*data(i,11));
end
data(:,19)=(data(:,16)+data(:,17))./data(:,18);

% output
cond=data(:,4)>10000;
figure
plot(data(cond,14),data(cond,15),'o')
hold
plot(data(cond,18),data(cond,19),'.')
hline(1)
xlabel("F7")
ylabel("(F3+F5)/F7)")
title(sprintf("%d Factor de linealidad (tras corregir con el deadtime)",brewer(brwid)))
legend(["op" "alt"])

figure
plot(data(cond,1),data(cond,15),'o')
hold
plot(data(cond,1),data(cond,19),'.')
hline(1)
xlabel("Date")
ylabel("(F3+F5)/F7)")
title(sprintf("%d Factor de linealidad (tras corregir con el deadtime)",brewer(brwid)))
legend(["op" "alt"])
% set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',1,'keeplimits','keepticks')




