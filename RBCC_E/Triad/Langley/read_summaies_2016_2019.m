
%path_root='D:\CODE\iberonesia\RBCC_E\Triad\2019\Triad'
%path_root='D:\CODE\iberonesia\RBCC_E\2019\Triad\langley'
%if ismac
% path_root='/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/Triad/2019/Triad'
%end

run(fullfile('..','read_config_'))

% remove year from cal.path_root
path_root=fullfile(Cal.path_root,'..');

ds_o=cell(3,1);
ds_a=ds_o;
ds=cell(3,2);
brewer=[157,183,185];
ano0=2014;

fblacklist=fullfile(path_root,'configs',['blacklist_',num2str(brewer(i)),'.txt']);  
confname=["Op","Alt"];
Width=20; Height=10; 



for i=1:3
    ds_o{i}=[];
    ds_a{i}=[];
    for ano=2015:2019
        j=ano-2014;
        %s1_=(strrep( strrep(fullfile(path_root,'summary_Brw157_2019.txt'),...
        %    '2019',num2str(ano)),'157',num2str(brewer(i))))
        s1_=strrep(fullfile(Cal.path_root,'Triad','Langley','triad_comp2019.mat'),'2019',num2str(ano));
        if exist(s1_)
          s=load(s1_,'A','ETC','F_corr','SL_B','SL_R','SL_corr_flag','summary','summary_old');
          t=write_summary_cfg(i,ano,s.summary_old,s.summary,s.SL_R,s.SL_B,s.A,s.ETC);
          t_sum(j,i,:)=t(i,:);
        else
          s=[];
        end
    end
end
 %  ds{i,1}=blacklist_summary(fblacklist,ds_o{i});
 ds_dep=cellfun(@(x) blacklist_summary(fblacklist,table2array(x)),t_sum,'UniformOutput',false);   
 %%Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6 A1 ETC
 label_= t_sum{1,1,1}.Properties.VariableNames
 
 for i=1:3
     ds{i,1}=cell2mat(ds_dep(:,i,1));
     ds{i,2}=cell2mat(ds_dep(:,i,2));
 end

 
% %%
% for j=1:2

j=1
 TSYNC=10; % (min)
 ref1=fix(ds{1,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
 ref2=fix(ds{2,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
 ref3=fix(ds{3,j}(:,1)*24*60/TSYNC)/24/60*TSYNC;
 s=scan_join([ref1,ds{1,j}],[ref2,ds{2,j}]);
 s=scan_join(s,[ref3,ds{3,j}]);
 ncol=18; % number of columns of the join matrix
 % % s 1 date
% %        2:19  157 'Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc A1 ETC idx'
% %       20:37  183
% %       38:55  185
% 
%% datos simultaneos ozono

 jsim=all(~isnan(s(:,7:ncol:end)),2) ; 
 ref=median(s(jsim,7:ncol:end),2);
 sim=s(jsim,:);
% 
%% prueba
mydata=sim(:,7:ncol:end)-mean(sim(:,7:ncol:end),2);
figure
plot(sim(:,1),mydata)
% %prueba
% 
label_= t_sum{1,1,1}.Properties.VariableNames

label_sim=['Date',cellfun(@(x) strcat(x,'_157'),label_,'UniformOutput',false),...
            cellfun(@(x) strcat(x,'_183'),label_,'UniformOutput',false),...
            cellfun(@(x) strcat(x,'_185'),label_,'UniformOutput',false)
            ];
 sim_table{j}=array2table(sortrows(s,1),'VariableNames',label_sim');
      
 
        %data=[ap_ev{i}.m(:,2:end),ap_ev{i}.std(:,2:end)];
        %l=size(data,2)/2;
        %idx=[1:l,;l+(1:l)];idx=idx(:);
        %data=data(:,idx);
        %lbl=[label_ap(7:end),cellfun(@(x) strcat(x,'_std'),label_ap(7:end),'UniformOutput',false)];
        %lbl=lbl(idx)
        %ap_ev_table{i}=array2table([ap_ev{i}.m(:,1),ap_ev{i}.N(:,2),data],'VariableNames',['Date','N',lbl])





%% medidas que difieren en menos de 1º
% histogram(sim(:,3:ncol:end)-mean(sim(:,3:ncol:end),2))
 jsim2=all(abs(sim(:,3:ncol:end)-mean(sim(:,3:ncol:end),2))<1,2);
 
 
 
 
%%
 ref1=median(sim(jsim2,7:ncol:end),2);  % reference 1st configuration  (O3_1)
 ref2=median(sim(jsim2,11:ncol:end),2);  % reference 2st configuration  (O3_2)
 ref3=median(sim(jsim2,13:ncol:end),2);  % reference 1st configuration+SL  (O3_1_sl)
% 
 ratio1=[sim(jsim2,1),100*(sim(jsim2,7:ncol:end)-ref1)./ref1];    % 1st configuration
 ratio2=[sim(jsim2,1),100*(sim(jsim2,11:ncol:end)-ref2)./ref2];   % 2nd  configurationn
 ratio3=[sim(jsim2,1),100*(sim(jsim2,13:ncol:end)-ref3)./ref3];   % 1st configuraiton + SL correction 
 
 m=[sim(jsim2,1),sim(jsim2,4:ncol:end)]; % airmas
 osc=m;
 osc(:,2:end)=m(:,2:end).*ref1 

 x1=smoothdata(ratio1(:,2:end),'gaussian',15,'SamplePoints',ratio1(:,1));
 x2=smoothdata(ratio2(:,2:end),'gaussian',15,'SamplePoints',ratio2(:,1));
 x3=smoothdata(ratio3(:,2:end),'gaussian',15,'SamplePoints',ratio3(:,1));
 
% f=figure;
% set(f,'Tag','comp_2016_2019');
 plot(ratio1(:,1),x1,'.'); hold all
% %plot(ratio2(:,1),x2,'-','LineWidth',2); hold all
% %plot(ratio2(:,1),x3,':'); hold all
 grid;
 datetick('x',12);
 title(['Configuracion',num2str(j)])
 title(strcat('Conf. ',confname(j)))
 set(gca,'Ylim',[-2,2]);
% end
% 
% ix=sort(findobj('-regexp','Tag','comp_2016_2019'));
% 
% printfiles_report(ix',Cal.dir_figs,'Width',Width,'Height',Height);
