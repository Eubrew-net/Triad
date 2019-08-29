% if ismac
%    path_root='/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/Triad/2019/Triad/ozone/';%Langley_Brw157_2019.txt'	
% end

read_config_;

path_root= '.'

%% eventos finales
ev{1}.dates=events{1}.dates([1,17,20,23]);
ev{1}.labels=[events{1}.labels([1,17,20,23])];

lgl=cell(3,1);
op=lgl;
alt=op;
lgl_ev=op;

brewer=[157,183,185];
ano0=2014;
for i=1:3
    lgl{i}=[];
    
        s1_=(strrep(fullfile(path_root,'Langley_Brw157_1619.txt'),...
                           '157',num2str(brewer(i))));
        if exist(s1_)
            s=load(s1_);
            lgl{i}=[lgl{i};s];
        else
           warning([ num2str(brewer(i)),'_',num2str(ano) ])
        end
        
    
    % unimos AM/PM
    lgl_o3{i}=sortrows([[lgl{i}(:,1)+0.25,lgl{i}(:,2:2:end)];[lgl{i}(:,1)+0.75,lgl{i}(:,3:2:end)]],1);
    op{i}=rows2vars(op_cfg{i});
    alt{i}=rows2vars(alt_cfg{i});
    
    %% por eventos
    data_tab_brw=meanperiods(lgl_o3{i}, events{i}); 
    data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1) ...
      round(data_tab_brw.m(:,4)) round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1) data_tab_brw.N(:,2)];
    lgl_ev{i}=array2table(data,'VariableNames',{'ETC1','err1','ETC2','err2','N'},'RowNames',str2name(data_tab_brw.evnts));
    
    %% ejemplo para tablas latex
%     input.data=lgl_ev{i};
%     % Switch transposing/pivoting your table if needed:
%     input.transposeTable = 1;
%     % Switch to generate a complete LaTex document or just a table:
%     input.makeCompleteLatexDocument = 0;
%     cellwrite('langley_eventos_.tex',latexTable(input));
%     
    
    %% operativa
    f=figure
    set(f,'Tag','ETC_S_OP');
    mean_smooth_abs(lgl_o3{i}(:,1),lgl_o3{i}(:,2),60,1);
    hold on
    [xx,yy]=stairs(op{i}.UsageDate,[op{i}.ETCOnO3Ratio]);
    h=plot([xx;[now]],[yy;yy(end,:)], '-','LineWidth',5);
    legend(h,'alt')
    %stairs(alt{i}.UsageDate,alt{i}.ETCOnO3Ratio,'k-','LineWidth',5)
    vline_v(events{i}.dates,' ',events{i}.labels)
    %title(num2str(brewer(i)));
    title([num2str(brewer(i)),' operativa']);
    grid on;
    datetick('x',12,'keeplimits','keepticks')
    xlim([datenum(2015,12,1),now])
    set(gca,'LooseInset',get(gca,'TightInset'))
    snapnow;
    
    
    %% Alternativa
    f=figure
    set(f,'Tag','ETC_S_ALT');
    mean_smooth_abs(lgl_o3{i}(:,1),lgl_o3{i}(:,4),60,1);
    hold on
    [xx,yy]=stairs(alt{i}.UsageDate,[alt{i}.ETCOnO3Ratio]);
    h=plot([xx;[now]],[yy;yy(end,:)], '-','LineWidth',5);
    legend(h,'op','alt')
    %stairs(alt{i}.UsageDate,alt{i}.ETCOnO3Ratio,'k-','LineWidth',5)
    vline_v(events{i}.dates,' ',events{i}.labels)
    title([num2str(brewer(i)),' alternativa']);
    grid on;
    datetick('x',12,'keeplimits','keepticks')
    xlim([datenum(2015,12,1),now])
    set(gca,'LooseInset',get(gca,'TightInset'))
    
end


% data_tab_brw=meanperiods(lgl_o3{i}, events{i}); 
% data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1) ...
%       round(data_tab_brw.m(:,4)) round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1) data_tab_brw.N(:,2)];
% lgl_ev_table{i}=array2table(data,'VariableNames',{'ETC1','err1','ETC2','err2','N'},'RowNames',str2name(data_tab_brw.evnts))

