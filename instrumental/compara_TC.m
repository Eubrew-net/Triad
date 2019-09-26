
s=load('IZO#185_sl_rw.mat')

config_temp.n_inst=1;
config_temp.brw_name="185";
config_temp.final_days=Cal.Date.FINAL_DAYS(1);

% first period
%config_temp.period=[datenum(2016,5,1),datenum(2017,12,31)];
%NTC=[0,-0.25,-0.17,-0.45,-0.57];

% Second period
config_temp.period=[datenum(2018,5,1),datenum(2019,3,15)];
NTC=[0,-0.23,-0.24,-0.18,-0.29];


[NTC1,ajuste1,Args1,Fraw,Fn1]=temp_coeff_raw(config_temp,s.sl_rw,'outlier_flag',1,'date_range',config_temp.period);
[NTC2,ajuste2,Args2,Fraw,Fn2]=temp_coeff_raw(config_temp,s.sl_rw,'outlier_flag',1,'N_TC',NTC','date_range',config_temp.period);



Forigx=Fn1; Fn=Fn2;
Cal.inst=1
Cal.brw_str{Cal.n_inst}=config_temp.brw_name;

% R6 -reference
figure;  set(gcf,'Tag','TEMP_COMP_DATE');
plot(Forigx(:,1),Forigx(:,2),'b.','MarkerSize',6); 
ylabel('Temperature','Color','b'); ax(1)=gca; set(ax(1),'YAxisLocation','right','XTicklabels',{' '}); 
[mn,sn]=grpstats(Forigx(:,[1,end,2]),{year(Forigx(:,1)),fix(Forigx(:,1))},{'mean','sem'});
[mt,st]=grpstats(Fn(:,[1,end]),{year(Fn(:,1)),fix(Fn(:,1))},{'mean','sem'}); 
ax(2) = axes('YAxisLocation','left','Color','none'); 
hold all; 
h1=errorbar(mn(:,1),mn(:,2),sn(:,2),'Color','k','Marker','s');
h2=errorbar(mt(:,1),mt(:,2),st(:,2),'Color','g','Marker','s');
title(['R6 Temperature dependence Brewer#', Cal.brw_str{Cal.n_inst}]); ylabel('Standard Lamp R6 ratio');
datetick('x',6,'KeepTicks','KeepLimits'); grid on; 
lg=legend([h1,h2],'Old temperature coeff','New temperature coeff','Location','best'); 
set(lg,'HandleVisibility','Off');  set(findobj(lg,'Type','text'),'FontSize',7,'HandleVisibility','Off');    
linkprop(ax,{'Position','XTick'}); 

figure; set(gcf,'Tag','TEMP_COMP_TEMP')
[mn,sn]=grpstats(Forigx(:,[2,end]),Forigx(:,2),{'mean','sem'});
[mt,st]=grpstats(Fn(:,[2,end]),Fn(:,2),{'mean','sem'});
hold on; 
h1=errorbar(mn(:,1),mn(:,2),sn(:,2),'Color','k','Marker','s','LineStyle','none');
plot(mn(:,1),mn(:,2),'Color','k','Marker','s','LineStyle','none')
h2=errorbar(mt(:,1),mt(:,2),st(:,2),'Color','g','Marker','s','LineStyle','none'); 
plot(mt(:,1),mt(:,2),'Color','g','Marker','s','LineStyle','none')
rline
grid;
lg=legend([h1,h2],{'Old temperature coeff','New temperature coeff'},'Location','best');
set(lg,'HandleVisibility','Off'); 
set(findobj(lg,'Type','text'),'FontSize',7,'HandleVisibility','Off');    
title(['R6 Temperature dependence Brewer#', Cal.brw_str{Cal.n_inst}]);
ylabel('Standard Lamp R6 ratio'); xlabel('Temperature'); set(gca,'Box','On');

boldify
