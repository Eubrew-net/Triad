clear all
run(fullfile('..','read_config_'))

% remove year from cal.path_root
path_root=fullfile(Cal.path_root,'..');

ds_o=cell(3,1);
ds_a=ds_o;
ds=cell(3,2);
brewer=[157,183,185];
ano0=2014;
flag={}; % depuration flag
fblacklist=fullfile(path_root,'configs',['blacklist_',num2str(brewer(i)),'.txt']);  
confname=["Op","Alt"];
Width=20; Height=10; 



for i=1:3
    ds_o{i}=[];
    ds_a{i}=[];
    for ano=2016:2019
        j=ano-2015;
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
    %%
    fblacklist=fullfile(path_root,'configs',['blacklist_',num2str(brewer(i)),'.txt'])
    ds_dep{i,:}=cellfun(@(x) blacklist_summary(fblacklist,table2array(x)),squeeze(t_sum(:,i,:)),'UniformOutput',false);
    
    
end
 % ds{i,1}=blacklist_summary(fblacklist,ds_o{i});
 %ds_dep=cellfun(@(x) blacklist_summary(fblacklist,table2array(x)),t_sum,'UniformOutput',false);   
 %%Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6 A1 ETC
 label_= t_sum{1,1,1}.Properties.VariableNames
 % bl=fullfile(Cal.path_root,'..','configs',['blacklist_',Cal.brw_str{1}],'.txt')
 % tbl=readtable(bl);tbl.Properties.VariableNames={'date_1','date_2','comment'};
 % x=ds_dep{:,:,1};y=t_sum{:,:,1};y.Date;j=ismember(y.Date,x(:,1));plot(y.Date,j);datetick
 ds_dep=cat(1,ds_dep);
 for i=1:3
     ds{i,1}=cell2mat(ds_dep{i}(:,1));
     ds{i,2}=cell2mat(ds_dep{i}(:,2));
     %raw{i}=table2array(vertcat(t_sum{:,i,1}));
     %x=ds{1,1}(:,1); y=raw{1}(:,1); j=ismember(y,x);plot(y,j,'-');
     %flag{i}=[y,j]; %depured data
 end

 
% %%
for j=1:2  % two configs

 TSYNC=10; % (min)
 ref1=fix(ds{1,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
 ref2=fix(ds{2,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
 ref3=fix(ds{3,j}(:,1)*24*60/TSYNC)/24/60*TSYNC;
 sj{j}=scan_join([ref1,ds{1,j}],[ref2,ds{2,j}]);
 sj{j}=scan_join(sj{j},[ref3,ds{3,j}]);
 ncol=18; % number of columns of the join matrix
 % % s 1 date
% %        2:19  157 'Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc A1 ETC idx'
% %       20:37  183
% %       38:55  185
% 
% 
label_= t_sum{1,1,1}.Properties.VariableNames
label_sim=['Date',cellfun(@(x) strcat(x,'_157'),label_,'UniformOutput',false),...
            cellfun(@(x) strcat(x,'_183'),label_,'UniformOutput',false),...
            cellfun(@(x) strcat(x,'_185'),label_,'UniformOutput',false)
            ];
 sim_table{j}=array2table(sortrows(sj{j},1),'VariableNames',label_sim');
 writetable(sim_table{1},strrep('sim_2016_2019_xx.csv','xx',confname{j}))



end

 
%% datos simultaneos ozono

 jsim=all(~isnan(sj{1}(:,7:ncol:end)),2) ; 
 ref=median(sj{1}(jsim,7:ncol:end),2);
 sim=sj{1}(jsim,:); 

%% tabla depuracion
   mydata=sim(:,7:ncol:end)-median(sim(:,7:ncol:end),2);
   
   figure
   plot(sim(:,1),mydata)
   datetick('x',12,'keeplimits')
   legend('157','183','185')
   
   
   % hourly means
   %[m,s]=grpstats([sim(:,1),mydata],fix(sim(:,1)*24)/24)
   % daily mean
   [m,sigma]=grpstats([sim(:,1),mydata],fix(sim(:,1)));

figure
subplot(2,1,1);
histogram(([sigma(:,2:end)]));
vline(0.4)
subplot(2,1,2);
histogram(([m(:,2:end)]));
vline(3)
sigma_out=abs(sigma(:,2:end))>0.4
media_out=abs(m(:,2:end))>3; %3

jx=any(sigma_out,2);
t_dep_day=array2table([m(jx,1),sigma_out(jx,:)]);
jx=any(media_out,2);
t_dep_day=vertcat(t_dep_day,array2table([m(jx,1),media_out(jx,:)]));


t_dep_day.Date=datetime(datestr(t_dep_day.Var1));
t_dep_day.Properties.VariableNames=['fecha',varname({Cal.brw_str{:}}),'Date'];
t_dep_day=table2timetable(t_dep_day);


writetable(timetable2table(t_dep_day))


 %%
        %data=[ap_ev{i}.m(:,2:end),ap_ev{i}.std(:,2:end)];
        %l=size(data,2)/2;
        %idx=[1:l,;l+(1:l)];idx=idx(:);
        %data=data(:,idx);
        %lbl=[label_ap(7:end),cellfun(@(x) strcat(x,'_std'),label_ap(7:end),'UniformOutput',false)];
        %lbl=lbl(idx)
        %ap_ev_table{i}=array2table([ap_ev{i}.m(:,1),ap_ev{i}.N(:,2),data],'VariableNames',['Date','N',lbl])





for j=1:2
 jsim=all(~isnan(sj{1}(:,7:ncol:end)),2) ; 
 ref=median(sj{1}(jsim,7:ncol:end),2);
 sim=sj{j}(jsim,:); 



% histogram(sim(:,3:ncol:end)-mean(sim(:,3:ncol:end),2))
% medidas que difieren en menos de 1º
jsim2=all(abs(sim(:,3:ncol:end)-mean(sim(:,3:ncol:end),2))<1,2);
 
%
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
 
 f=figure;
% set(f,'Tag','comp_2016_2019');
 plot(ratio1(:,1),x1,'.'); hold all
% %plot(ratio2(:,1),x2,'-','LineWidth',2); hold all
% %plot(ratio2(:,1),x3,':'); hold all
 grid;
 datetick('x',12);
 title(['Configuracion',num2str(j)])
 title(strcat('Conf. ',confname(j)))
 set(gca,'Ylim',[-2,2]);
 legend('157','183','185')
 
end
% end
% 
% ix=sort(findobj('-regexp','Tag','comp_2016_2019'));
% 
% printfiles_report(ix',Cal.dir_figs,'Width',Width,'Height',Height);
