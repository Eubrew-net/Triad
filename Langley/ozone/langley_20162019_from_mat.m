clear all
close all
% ejecutamos la configuracion
%run(fullfile(path_root,'..','Triad','read_config_'));
run(fullfile('..','..','read_config_'));
path_root=fullfile(Cal.path_root);
% not used
%fblacklist=fullfile(path_root,'configs',['blacklist_',num2str(brewer(i)),'.txt']); 


%% eventos langley finales
% 157
ev=cell(3,1);
%%
lgl=cell(3,1);
op=lgl;
alt=op;
lgl_ev=op;

brewer=[157,183,185];
ano0=2014;
s=[];


%% eventos langley finales
ev=cell(3,1);

% 157
importantes{1}=[1,4,16,20,23];
ev{1}.dates=events{1}.dates(importantes{1})';
ev{1}.labels=events{1}.labels(importantes{1});
% a�adimos ios 2015
ev{1}.dates(1)=datenum(2015,6,6);
ev{1}.labels(1)={'IOS'};
ev{1}.op_ETC=[op_cfg{1}{8,[importantes{1}]}];
ev{1}.alt_ETC=[alt_cfg{1}{8,importantes{1}}];
ev{1}.op_A1=[op_cfg{1}{7,importantes{1}}];
ev{1}.alt_A1=[alt_cfg{1}{7,importantes{1}}];
ev{1}.fechas=datetime(datestr(ev{1}.dates))

%% remove 

% 183
importantes{2}=[1,9,13,15,16];
ev{2}.dates=events{2}.dates(importantes{2})';
ev{2}.labels=events{2}.labels(importantes{2});
% a�adimos ios 2015
ev{2}.dates(1)=datenum(2015,6,9);
ev{2}.labels(1)={'IOS'};
ev{2}.op_ETC=[op_cfg{2}{8,[importantes{2}]}];
ev{2}.alt_ETC=[alt_cfg{2}{8,importantes{2}}];
ev{2}.op_A2=[op_cfg{2}{7,importantes{2}}];
ev{2}.alt_A2=[alt_cfg{2}{7,importantes{2}}];
ev{2}.fechas=datetime(datestr(ev{2}.dates))';
% 185
%importantes{3}=[4,7,8,10,14,16,18,19,24,26,28,30,31];
importantes{3}=[1,4,7,8,12,17,18,19,20,26,27,29,31,32];
ev{3}.dates=events{3}.dates(importantes{3})';
ev{3}.labels=events{3}.labels(importantes{3});
% a�adimos ios 2015
ev{3}.dates(1)=datenum(2015,6,10);
ev{3}.labels(1)={'IOS'};
ev{3}.op_ETC=[op_cfg{3}{8,[importantes{3}]}];
ev{3}.alt_ETC=[alt_cfg{3}{8,importantes{3}}];
ev{3}.op_A3=[op_cfg{3}{7,importantes{3}}];
ev{3}.alt_A3=[alt_cfg{3}{7,importantes{3}}]
ev{3}.fechas=datetime(datestr(ev{3}.dates))';

table(ev{3}.fechas',ev{3}.labels')




%% load data

for i=1:3
    
    %%
    lgl{i}=[];
    
     s1_=(strrep('Langley_Brw157_1619.txt','157',num2str(brewer(i))));
     if exist(s1_)
            s=load(s1_);
            lgl{i}=[lgl{i};s];
     else          
           disp(s1_)
     end
        

 
%%
    % unimos AM/PM
    %Date etc_1_am etc_1_pm slope_1_am slope_1_pm  etc_2_am etc_2_pm slope_2_am slope_2_pm detc_1_am detc_1_pm dslope_1_am dslope_1_pm  detc_2_am detc_2_pm dslope_2_am dslope_2_p
    lgl_o3{i}=sortrows([[lgl{i}(:,1)+0.25,lgl{i}(:,2:2:end)];[lgl{i}(:,1)+0.75,lgl{i}(:,3:2:end)]],1);
    %Date etc_1 slope_1  etc_2 slope_2 detc_1 dslope_1  detc_2 dslope_2
    op{i}=rows2vars(op_cfg{i});
    alt{i}=rows2vars(alt_cfg{i});

% end
    %% por eventos
    data_tab_brw=meanperiods(lgl_o3{i}, events{i}); 
    data=[round(data_tab_brw.m(:,2)), round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1), ...
          op_cfg{i}{8,:}',...
      round(data_tab_brw.m(:,4)), round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1),...
          alt_cfg{i}{8,:}',...
      round(data_tab_brw.m(:,[6])), round(data_tab_brw.std(:,[6])./sqrt(data_tab_brw.N(:,2)),1),...
      round(data_tab_brw.m(:,[8])), round(data_tab_brw.std(:,[8])./sqrt(data_tab_brw.N(:,2)),1),...
      data_tab_brw.N(:,2)
               
          ];
    
    lgl_ev{i}=array2table(data,'VariableNames',{'ETC1','err1','ETC_op','ETC2','err2','ETC_alt','dETC1','derr1','dETC2','derr2','N'},'RowNames',str2name(data_tab_brw.evnts));
    lgl_ev{i}.Date=events{i}.dates;
       
    %lgl_ev{i}.Date=datetime(datestr(events{i}.dates));
    %lgl_ev{i}=timetable2table(table2timetable(lgl_ev{i}));
    
    %% eventos finales
    if ~isempty(ev{i})
       data_tab_brw=meanperiods(lgl_o3{i}, ev{i}); 
       %data=[round(data_tab_brw.m(:,2)), round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1), ev{i}.op_ETC', ...
       %round(data_tab_brw.m(:,4)), round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1), ev{i}.alt_ETC', data_tab_brw.N(:,2)];
       data=[round(data_tab_brw.m(:,2)), round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1), ...
          ev{i}.op_ETC',...
      round(data_tab_brw.m(:,4)), round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1),...
           ev{i}.alt_ETC',...
      round(data_tab_brw.m(:,[6])), round(data_tab_brw.std(:,[6])./sqrt(data_tab_brw.N(:,2)),1),...
      round(data_tab_brw.m(:,[8])), round(data_tab_brw.std(:,[8])./sqrt(data_tab_brw.N(:,2)),1),...
      data_tab_brw.N(:,2)
               
          ];
       
       
       lgl_ev_final{i}=array2table(data,'VariableNames',{'ETC1','err1','ETC_op','ETC2','err2','ETC_alt','dETC1','derr1','dETC2','derr2','N'},'RowNames',str2name(data_tab_brw.evnts));
       %array2table(data,'VariableNames',{'ETC1','err1','ETC_op','ETC2','err2','ETC_alt','N'},'RowNames',str2name(data_tab_brw.evnts));
       %lgl_ev_final{i}.Date=str2num(datestr(ev{i}.dates,'yyyymmdd'))
       lgl_ev_final{i}.Date=datetime(datestr(ev{i}.dates));
       lgl_ev_final{i}=timetable2table(table2timetable(lgl_ev_final{i}));
       lgl_ev_final{i}=addvars(lgl_ev_final{i},alt_cfg{i}{7,importantes{i}}','After','N','NewVariableName','A1_chk');
       lgl_ev_final{i}=addvars(lgl_ev_final{i},op_cfg{i}{7,importantes{i}}','After','N','NewVariableName','A1_op');
       
       
       disp(lgl_ev_final{i})
       writetable(lgl_ev_final{i},'IzoTriad_2016_2019.xls','Sheet',strcat('ETC_final',Cal.brw_str{i}))
    end
    
    % ejemplo para tablas latex
    input.data=lgl_ev{i};
    % Switch transposing/pivoting your table if needed:
    input.transposeTable = 0;
    % Switch to generate a complete LaTex document or just a table:
    %input.makeCompleteLatexDocument = 0;
    %l=latexTable(input);
    %cellwrite(sprintf('langley_eventos_%03d.tex',brewer(i)),l);
    %%
    
    %% outputs 
    figure 
    % lgl_o3
    %Date etc_1 slope_1  etc_2 slope_2 detc_1 dslope_1  detc_2 dslope_2
  
    [h1,p1]=mean_smooth_abs(lgl_o3{i}(:,1),lgl_o3{i}(:,6),60,0);hold on
    [h2,p2]=mean_smooth_abs(lgl_o3{i}(:,1),lgl_o3{i}(:,8),60,0);hold on
    
    px=[p1(:,1:2),p1(:,2)+[-p1(:,3),p1(:,3)]];h1=ploty(px,'-')
    px=[p2(:,1:2),p2(:,2)+[-p2(:,3),p2(:,3)]];h2=ploty(px,'-.')
    hold on
    plot(lgl_o3{i}(:,1),lgl_o3{i}(:,6),'.')
    plot(lgl_o3{i}(:,1),lgl_o3{i}(:,8),'x')
    
    ylim([prctile(lgl_o3{i}(:,6),0.1),prctile(lgl_o3{i}(:,6),99.9)])
    
    [xx,yy]=stairs(datenum(lgl_ev_final{i}.Date),[lgl_ev_final{i}.dETC1,lgl_ev_final{i}.dETC2]);
    h=plot([xx;[now,now]],[yy;yy(end,:)], '-','LineWidth',5);
    legend(h,'op','alt')
    %stairs(alt{i}.UsageDate,alt{i}.ETCOnO3Ratio,'k-','LineWidth',5)
    %vline_v(events{i}.dates,' ',events{i}.labels)
    vline_v(ev{i}.dates,' ',str2latex(ev{i}.labels))
    title(num2str(brewer(i)));
    %title([num2str(brewer(i)),'operativa']);
    grid on;
    datetick('x',12,'keeplimits','keepticks')
    xlim([datenum(2015,12,1),now])
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    %% outputs 
    figure 
    mean_smooth_abs(lgl_o3{i}(:,1),lgl_o3{i}(:,2),60,1);
    hold on
    plot(lgl_o3{i}(:,1),lgl_o3{i}(:,2),'.')
    ylim([prctile(lgl_o3{i}(:,2),0.1),prctile(lgl_o3{i}(:,2),99.9)])
    [xx,yy]=stairs(op{i}.UsageDate,[op{i}.ETCOnO3Ratio,alt{i}.ETCOnO3Ratio]);
    h=plot([xx;[now,now]],[yy;yy(end,:)], '-','LineWidth',5);
    legend(h,'op','alt')
    %stairs(alt{i}.UsageDate,alt{i}.ETCOnO3Ratio,'k-','LineWidth',5)
    %vline_v(events{i}.dates,' ',events{i}.labels)
    vline_v(ev{i}.dates,' ',str2latex(ev{i}.labels))
    %title(num2str(brewer(i)));
    title([num2str(brewer(i)),'operativa']);
    grid on;
    datetick('x',12,'keeplimits','keepticks')
    %xlim([datenum(2015,12,1),now])
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    
    %% operative
    figure
    
    mean_smooth_abs(lgl_o3{i}(:,1),lgl_o3{i}(:,4),60,1);
    hold on
    plot(lgl_o3{i}(:,1),lgl_o3{i}(:,4),'.')
    ylim([prctile(lgl_o3{i}(:,4),0.1),prctile(lgl_o3{i}(:,4),99.9)])
    [xx,yy]=stairs(op{i}.UsageDate,[op{i}.ETCOnO3Ratio,alt{i}.ETCOnO3Ratio]);
    h=plot([xx;[now,now]],[yy;yy(end,:)], '-','LineWidth',5);
    legend(h,'op','alt')
    %stairs(alt{i}.UsageDate,alt{i}.ETCOnO3Ratio,'k-','LineWidth',5)
    %vline_v(events{i}.dates,' ',events{i}.labels)
 	vline_v(ev{i}.dates,' ',str2latex(ev{i}.labels))
    title([num2str(brewer(i)),' alternativa']);
    grid on;
    datetick('x',12,'keeplimits','keepticks')
    %xlim([datenum(2015,12,1),now])
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    %% summary
    disp(lgl_ev_final{i})
    
    snapnow
 end


% data_tab_brw=meanperiods(lgl_o3{i}, events{i}); 
% data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1) ...
%       round(data_tab_brw.m(:,4)) round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1) data_tab_brw.N(:,2)];
% lgl_ev_table{i}=array2table(data,'VariableNames',{'ETC1','err1','ETC2','err2','N'},'RowNames',str2name(data_tab_brw.evnts))
