
clear all;
run(fullfile('..','read_config_'));

brewer=Cal.brw
for i=3:-1:1
    
    Cal.n_inst=i
    tablasc{i}=report_sc(Cal,Cal.brw_config_files{i,1},'grp','events');
    %tabla_sc{i}=display_table(tablasc{i}.data(:,2:end),tablasc{i}.data_lbl,15,'.0f',tablasc{i}.events);
    
    
    %jn=~isnan(tablasc{i}.data(:,2)); % para quitar eventos sin datos
    jn=all(tablasc{i}.data(:,1),2);
    tabla_sc{i}=array2table(tablasc{i}.data(jn,1:end),'VariableNames',['Date',varname(tablasc{i}.data_lbl)],'Rownames',strrep(varname(tablasc{i}.events(jn)),'x_',''));
    tabla_sc{i}.Date=datetime(datestr(tablasc{i}.data(jn,1)));
    %writetable(dsp_ev{i},strrep('dsp_ev_157_2016_2019.txt','157',num2str(brewer(i))))
    tabla_sc{i}=table2timetable(tabla_sc{i});
    disp(tabla_sc{i})
    tabla_sc{i}=timetable2table(tabla_sc{i});
    writetable(tabla_sc{i},'IzoTriad_2016_2019.xls','Sheet',strrep('sc_avg_157','157',num2str(brewer(i))),'WriteRowNames',true)
    
    
    %writetable(tabla_sc_157,'SC_summary_157.txt')
    %writetable(tabla_sc_157,'SC_summary_157.xls')
    
    %% figura
    figure; set(gcf,'Tag',strcat('SC_CSN_hist',Cal.brw_str{Cal.n_inst}));
    patch([min(tablasc{i}.data(:,1))-10 max(tablasc{i}.data(:,1))+10 max(tablasc{i}.data(:,1))+10 min(tablasc{i}.data(:,1))-10],...
        [repmat(1025,1,2) repmat(1027,1,2)], ...
        [.953,.953,.953],'LineStyle',':'); hold on;
    errorbar(tablasc{i}.data(:,1),tablasc{i}.data(:,2),matadd(tablasc{i}.data(:,2),-tablasc{i}.data(:,3)),matadd(tablasc{i}.data(:,4),-tablasc{i}.data(:,2)),'*');
    set(gca,'Layer','Top');
    grid; datetick('x',12,'keepLimits','keepTicks');
    title(sprintf('Brewer %s ',Cal.brw_name{Cal.n_inst}));
    ylabel('Cal Step Number'); box on; axis('tight');
    hold on; s=stairs(tablasc{i}.data(:,1),tablasc{i}.data(:,end-1),'-','color','m','LineWidth',2);
    vl=vline_v(op_cfg{i}{1,:}','-r', op_cfg{i}.Properties.VariableNames');
    %hline=
    set(vl,'LineStyle','None'); set(findobj(gcf,'Type','Text'),'FontSize',7);
    
end

