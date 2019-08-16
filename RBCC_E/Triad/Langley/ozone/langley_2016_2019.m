%read_config_;

lgl=cell(3,1);
op=lgl;
alt=op;
lgl_ev=op;

brewer=[157,183,185];
ano0=2014;
for i=1:3
    lgl{i}=[];
    
    for ano=2015:2019
        j=ano-2014;
        s1_=(strrep( strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley/ozone/Langley_Brw157_2019.txt','2019',num2str(ano)),'157',num2str(brewer(i))))
        if exist(s1_)
            s=load(s1_);
            lgl{i}=[lgl{i};s];
        end
        
    end
    % unimos AM/PM
    lgl_o3{i}=sortrows([[lgl{i}(:,1)+0.25,lgl{i}(:,2:2:end)];[lgl{i}(:,1)+0.75,lgl{i}(:,3:2:end)]],1)
    op{i}=rows2vars(op_cfg{i});
    alt{i}=rows2vars(alt_cfg{i});
    
    %% por eventos
    data_tab_brw=meanperiods(lgl_o3{i}, events{i}); 
    data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1) ...
      round(data_tab_brw.m(:,4)) round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1) data_tab_brw.N(:,2)];
    lgl_ev{i}=array2table(data,'VariableNames',{'ETC1','err1','ETC2','err2','N'},'RowNames',str2name(data_tab_brw.evnts));

    
    
    figure(i)
    
    mean_smooth_abs(lgl_o3{i}(:,1),lgl_o3{i}(:,2),60,1)
    hold on
    [xx,yy]=stairs(op{i}.UsageDate,[op{i}.ETCOnO3Ratio,alt{i}.ETCOnO3Ratio]);
    plot([xx;[now,now]],[yy;yy(end,:)], '-','LineWidth',5);
    %stairs(alt{i}.UsageDate,alt{i}.ETCOnO3Ratio,'k-','LineWidth',5)
    vline_v(events{i}.dates,' ',events{i}.labels)
    title(num2str(brewer(i)));
    grid on;
    datetick('keepticks')
    
end

% data_tab_brw=meanperiods(lgl_o3{i}, events{i}); 
% data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1) ...
%       round(data_tab_brw.m(:,4)) round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1) data_tab_brw.N(:,2)];
% lgl_ev_table{i}=array2table(data,'VariableNames',{'ETC1','err1','ETC2','err2','N'},'RowNames',str2name(data_tab_brw.evnts))

