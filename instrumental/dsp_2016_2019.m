% collecting dsp information
clear all
read_config_
%Cal.dir_figs=fullfile('..',Cal.dir_figs)
%Cal.dir_tables=fullfile('..',Cal.dir_tables)
brewer=Cal.brw

% Tabla por eventos
tabla_dsp_ev=cell(3,1);

% tabla general
tabla_dsp=tabla_dsp_ev;

% sumario dispersion cuadratico
dsp_quad=cell(3,1);
% sumario dispersion c�bico
dsp_cubic=cell(3,1);

% resumen por eventos
dsp_ev=cell(3,1);
% resumen general
dsp_table=cell(3,1);


for i=1:3

% salidas del report    
tabla_dsp_ev{i}=[]; tabla_dsp{i}=[];

% sumarios
dsp_cubic{i}=[];dsp_cuad{i}=[];

% resumenes
dsp_ev{i}=[];dsp_table{i}=[];

Cal.Date.CALC_DAYS=datenum(2015,6,1):now;
Cal.n_inst=i
%average for events
% same events as config
tabla_dsp_ev{i}=report_dispersion_new(Cal,'grp_custom',events{i},'fpath',fullfile(Cal.path_root,'..','DSP'),'process',0);%,...
                                     %      'date_range',Cal.Date.CALC_DAYS);  
                                       
% weekly average in practice all data                                       
[tabla_dsp{i},dsp_quad{i},dsp_cubic{i}]=report_dispersion_new(Cal,'grp','week','fpath',fullfile(Cal.path_root,'..','DSP'),'process',0,...
                                           'date_range',Cal.Date.CALC_DAYS);  
%% tabla sumario 
jn=~isnan(tabla_dsp{i}.data(:,2));
dsp_table{i}=array2table(tabla_dsp{i}.data(jn,1:end),'VariableNames',varname(tabla_dsp{i}.data_lbl));
dsp_table{i}.Date=datetime(datestr(tabla_dsp{i}.data(jn,1)));
%writetable(dsp_table{i},strrep('dsp_157_2016_2019.txt','157',num2str(brewer(i))))
dsp_table{i}=table2timetable(dsp_table{i});
disp(dsp_table{i})
dsp_table{i}=timetable2table(dsp_table{i});
writetable(dsp_table{i},'IzoTriad_2016_2019.xls','Sheet',strrep('dsp_157','157',num2str(brewer(i))),'WriteRowNames',true)
%strrep('dsp_157_2016_2019.txt','157',num2str(brewer(i))))
%writetable(dt_evt{i},'IzoTriad_2016_2019.xls','Sheet',strrep('dt_157','157',num2str(brewer(i))),'WriteRowNames',true)
%% tabla por eventos 
%A1_op=op_cfg{i}(7,:);
%strmatch('std',tabla_dsp_ev{i}.data_lbl)
% tabla_dsp_ev{i}.data(:,16)=op_cfg{i}{7,:};
% tabla_dsp_ev{i}.data_lbl{16}='A1_op';
% tabla_dsp_ev{i}.data(:,18)=alt_cfg{i}{7,:};
% tabla_dsp_ev{i}.data_lbl{18}='A1_chk';
%jn=~isnan(tabla_dsp_ev{i}.data(:,2));
jn=all(tabla_dsp_ev{i}.data(:,1),2);
dsp_ev{i}=array2table(tabla_dsp_ev{i}.data(jn,1:end),'VariableNames',varname(tabla_dsp_ev{i}.data_lbl),'Rownames',strrep(varname(tabla_dsp_ev{i}.events(jn)),'x_',''));
dsp_ev{i}.Date=datetime(datestr(tabla_dsp_ev{i}.data(jn,1)));
dsp_ev{i}=addvars(dsp_ev{i},alt_cfg{i}{7,:}','After','CSN','NewVariableName','A1_chk');
dsp_ev{i}=addvars(dsp_ev{i},op_cfg{i}{7,:}','After','CSN','NewVariableName','A1_op');
dsp_ev{i}=movevars(dsp_ev{i},{'A1Quad_','A1Cubic'},'After','CSN');

%writetable(dsp_ev{i},strrep('dsp_ev_157_2016_2019.txt','157',num2str(brewer(i))))
dsp_ev{i}=table2timetable(dsp_ev{i});
jn=~isnan(tabla_dsp_ev{i}.data(:,2));
disp(dsp_ev{i}(jn,:))
dsp_ev{i}=timetable2table(dsp_ev{i}(jn,:));

writetable(dsp_ev{i},'IzoTriad_2016_2019.xls','Sheet',strrep('dsp_avg_157','157',num2str(brewer(i))),'WriteRowNames',true)
%%
f1=figure
 dx=dsp_ev{i};
h1=plot(dx.Date,dx.A1Quad_-dx.A1_op,'-rx',dx.Date,dx.A1Cubic-dx.A1_op,'-ro')
hold all
h1=plot(dx.Date,dx.A1Quad_-dx.A1_chk,':bx',dx.Date,dx.A1Cubic-dx.A1_chk,':bo')
grid
title(Cal.brw_str{i})



%% figura DSP A1 Residuos frente al promedio -> deberia ser frente a la configuracion
f1=figure
set(f1,'Tag',strcat('DSP_A1_RES',Cal.brw_str{i})); 

mq{i}=mean_smooth_abs(tabla_dsp{i}.data(:,1),tabla_dsp{i}.data(:,[15]),90,0);
mc{i}=mean_smooth_abs(tabla_dsp{i}.data(:,1),tabla_dsp{i}.data(:,[17]),90,0);
hold on
m{i}=nanmean(tabla_dsp{i}.data(:,[15,17]))
s{i}=nanstd(tabla_dsp{i}.data(:,[15,17]))
n{i}=sum(~isnan(tabla_dsp{i}.data(:,[15,17])));
e{i}=nanstd(tabla_dsp{i}.data(:,[15,17]))./sqrt(n{i})
h1=plot(tabla_dsp{i}.data(:,1),(mq{i}(:,1)-m{i}(1))/0.001)
hold on
h2=plot(tabla_dsp{i}.data(:,1),(mc{i}(:,1)-m{i}(2))/0.001)

plot(tabla_dsp{i}.data(:,1),(tabla_dsp{i}.data(:,15)-m{i}(1))/0.001,'k+')
hold on
plot(tabla_dsp{i}.data(:,1),(tabla_dsp{i}.data(:,17)-m{i}(2))/0.001,'bx')

datetick('x','mmm/yy','keepticks')
grid on
ylabel('Micromenter step')
title({strrep('Brewer #157 Ozone Abs coeff diff vs mean','157',num2str(brewer(i))),sprintf('quad=%.4f cubic=%.4f',m{i})})
box on
set(gca,'YLim',[-2 2])
vline_v(events{i}.dates,'k',events{i}.labels);
hfill([-0.5,0.5],'grey')
legend([h1,h2],{'cuadratic','cubic'})



%% figura DSP hist�rica
                                

f2=figure;
set(f2,'Tag',strcat('DSP_A1_HIST',Cal.brw_str{i})); 
tabla_dsp{i}.data(isnan(tabla_dsp{i}.data(:,2)),:)=[];
l1=ploty(tabla_dsp{i}.data(:,[1 15 17]),'s');
set(l1(1),'Marker','s','LineStyle','-','Color','k','MarkerFaceColor','k'); 
set(l1(2),'Marker','s','LineStyle','-','Color','g','MarkerFaceColor','g');
%lim=get(gca,'XLim'); set(gca,'YLim',[0.3360 0.3435],'XLim',[lim(1)-10 lim(2)+10]); 
set(gca,'YLim',round(m{i},4)+[-0.005,0.005])
hold on; st(2)=stairs([events_cfg_chk{i}.data(1,:) tabla_dsp{i}.data(end,1)],[events_cfg_chk{i}.data(7,:) events_cfg_chk{i}.data(7,end)],'-','color','m','LineWidth',2);
         st(1)=stairs([events_cfg_op{i}.data(1,:) tabla_dsp{i}.data(end,1)],[events_cfg_op{i}.data(7,:) events_cfg_op{i}.data(7,end)],'-','color','b','LineWidth',2);
vl=vline_v(events_cfg_op{i}.data(1,:),'-r',Cal.events_raw{Cal.n_inst}(cell2mat(Cal.events_raw{Cal.n_inst}(:,2))>=events_cfg_op{i}.data(1,1),3)); 
set(vl,'LineStyle','None'); 
set(findobj(gca,'Type','text'),'FontSize',7);    
title(sprintf('A1: %s',Cal.brw_name{i}));
datetick('x',12,'KeepLimits','keepTicks'); ylabel('Ozone Absorption Coefficient (A1)'); grid
legendflex([l1;st(1);st(2)],{'A1 Quad','A1 Cubic','A1 ref (Op)','A1 ref (Chk)'}, ...
                           'anchor', {'sw','sw'}, 'buffer',[10 -5], 'nrow',2,'fontsize',8,'box','on');
set(findobj(gcf,'Type','text'),'BackgroundColor','w');

printfiles_report(gcf,Cal.dir_figs,'Width',21,'Height',9,'LineMode','Scaled');



end