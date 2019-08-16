%read_config_;

ds_o=cell(3,1);
ds_a=ds_o;


brewer=[157,183,185];
ano0=2014;
for i=1:3
    ds_o{i}=[];
    ds_a{i}=[];
    for ano=2015:2019
        j=ano-2014;
        s1_=(strrep( strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley/summary_Brw157_2019.txt','2019',num2str(ano)),'157',num2str(brewer(i))))
        s2_=(strrep( strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley/summary_old_Brw157_2019.txt','2019',num2str(ano)),'157',num2str(brewer(i))))
       
        if exist(s1_)
            s=load(s1_);
            ds_o{i}=[ds_o{i};s(:,1:15)];
        end
        
        if exist(s2_)
            s=load(s2_);
            ds_a{i}=[ds_a{i};s(:,1:15)];
        end
    end
    
    %%Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc
    
    
%     % unimos AM/PM
%     lgl_o3{i}=sortrows([[lgl{i}(:,1)+0.25,lgl{i}(:,2:2:end)];[lgl{i}(:,1)+0.75,lgl{i}(:,3:2:end)]],1)
%     op{i}=rows2vars(op_cfg{i});
%     alt{i}=rows2vars(alt_cfg{i});
%     
%     %% por eventos
%     data_tab_brw=meanperiods(lgl_o3{i}, events{i}); 
%     data=[round(data_tab_brw.m(:,2)) round(data_tab_brw.std(:,2)./sqrt(data_tab_brw.N(:,2)),1) ...
%       round(data_tab_brw.m(:,4)) round(data_tab_brw.std(:,4)./sqrt(data_tab_brw.N(:,2)),1) data_tab_brw.N(:,2)];
%     lgl_ev{i}=array2table(data,'VariableNames',{'ETC1','err1','ETC2','err2','N'},'RowNames',str2name(data_tab_brw.evnts));

    
    
%     figure(i)
%     
%     mean_smooth_abs(lgl_o3{i}(:,1),lgl_o3{i}(:,2),60,1)
%     hold on
%     [xx,yy]=stairs(op{i}.UsageDate,[op{i}.ETCOnO3Ratio,alt{i}.ETCOnO3Ratio]);
%     plot([xx;[now,now]],[yy;yy(end,:)], '-','LineWidth',5);
%     %stairs(alt{i}.UsageDate,alt{i}.ETCOnO3Ratio,'k-','LineWidth',5)
%     vline_v(events{i}.dates,' ',events{i}.labels)
%     title(num2str(brewer(i)));
%     grid on;
%     datetick('keepticks')
    
end

TSYNC=10; % (min)
ref1=fix(ds_o{1}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
ref2=fix(ds_o{2}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
ref3=fix(ds_o{3}(:,1)*24*60/TSYNC)/24/60*TSYNC;
s=scan_join([ref1,ds_o{1}],[ref2,ds_o{2}]);
s=scan_join(s,[ref3,ds_o{3}]);
%%Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc
jsim=all(~isnan(s(:,7:15:end)),2) ; 
ref=median(s(jsim,7:15:end),2);
sim=s(jsim,:);
% medidas que difieren en menos de 1º
%histogram(sim(:,3:15:end)-mean(sim(:,3:15:end),2))
jsim2=all(abs(sim(:,3:15:end)-mean(sim(:,3:15:end),2))<1.0,2);
ref=median(sim(jsim2,7:15:end),2);

ratio=[sim(jsim2,1),100*(sim(jsim2,7:15:end)-ref)./ref];
ratio2=[sim(jsim2,1),100*(sim(jsim2,11:15:end)-ref)./ref];
ratio3=[sim(jsim2,1),100*(sim(jsim2,13:15:end)-ref)./ref];

m=[sim(jsim2,1),sim(jsim2,4:15:end)];

x=smoothdata(ratio(:,2:end),'gaussian',15,'SamplePoints',ratio(:,1));
x2=smoothdata(ratio2(:,2:end),'gaussian',15,'SamplePoints',ratio(:,1));
x3=smoothdata(ratio3(:,2:end),'gaussian',15,'SamplePoints',ratio(:,1));
figure
plot(ratio(:,1),x,'.'); hold all
%plot(ratio(:,1),x2,'-'); hold all
plot(ratio(:,1),x3,':'); hold all
grid;
