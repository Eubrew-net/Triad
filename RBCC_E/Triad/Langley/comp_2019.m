


clear all;  close all;
run(fullfile('..','read_config_'));


% remove year from cal.path_root
%path_root=Cal.path_root(1:(end-5));
path_root=fullfile(Cal.path_root,'..');

confname=["Op","Alt"];
 
ds_o=cell(3,1);
ds_a=ds_o;

brewer=[157,183,185];
%% eventos 2019
bevents=brw_events(year(brw_events.Date)==2019,:)

for i=1:3
    ds_o{i}=[];
    ds_a{i}=[];
    for ano=2019:2019
        j=ano-2018;
        %s1_=(strrep( strrep(fullfile(path_root,'summary_Brw157_2019.txt'),...
        %    '2019',num2str(ano)),'157',num2str(brewer(i))))
        s1_=fullfile(path_root,'Triad',num2str(ano),'Triad',['summary_Brw',num2str(brewer(i)),'_',num2str(ano),'.txt']);
        %s2_=(strrep( strrep(fullfile(path_root,'summary_old_Brw157_2019.txt'),...
        %    '2019',num2str(ano)),'157',num2str(brewer(i))));
        s2_=fullfile(path_root,'Triad',num2str(ano),'Triad',['summary_old_Brw',num2str(brewer(i)),'_',num2str(ano),'.txt']);
       
        if exist(s1_)
            s=load(s1_);
            ds_a{i}=[ds_a{i};s(:,1:15)];
        else
           warning([s1_ ,'not found'])
        end    
        
        if exist(s2_)
            s=load(s2_);
            ds_o{i}=[ds_o{i};s(:,1:15)];
        else
         warning([s2_ ,'not found'])
        end
        
        fprintf('Año %d brewer %d diff(o1,a2)=%d diff(o1,a2)=%d\n',ano,brewer(i),sum((ds_o{i}(:,6)-ds_a{i}(:,10))~=0),sum((ds_o{i}(:,10)-ds_a{i}(:,6))~=0));
    end
    
    %%Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc
    
 fblacklist=fullfile(path_root,'configs',['blacklist_',num2str(brewer(i)),'.txt']);  
    
 ds{i,1}=blacklist_summary(fblacklist,ds_o{i});
 ds{i,2}=blacklist_summary(fblacklist,ds_a{i});
     
end

%% Comparamos los 3 instrumentos
for j=1:2
 
TSYNC=10; % (min)
ref1=fix(ds{1,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
ref2=fix(ds{2,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
ref3=fix(ds{3,j}(:,1)*24*60/TSYNC)/24/60*TSYNC;
s=scan_join([ref1,ds{1,j}],[ref2,ds{2,j}]);
s=scan_join(s,[ref3,ds{3,j}]);
label{1}='Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc';
% s 1 date
%        2:16  157 'Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc'
%       17:31  183
%       32:45  185

% datos simultaneos
jsim=all(~isnan(s(:,7:15:end)),2) ; 
ref=median(s(jsim,7:15:end),2);
sim=s(jsim,:);

%histogram(sim(:,3:15:end)-mean(sim(:,3:15:end),2))
jsim2=all(abs(sim(:,3:15:end)-mean(sim(:,3:15:end),2))<1.0,2);

ref1=median(sim(jsim2,7:15:end),2);  % reference 1st configuration  (O3_1)
ref2=median(sim(jsim2,11:15:end),2);  % reference 2st configuration  (O3_2)
ref3=median(sim(jsim2,13:15:end),2);  % reference 1st configuration+SL  (O3_1_sl)

ratio1=[sim(jsim2,1),100*(sim(jsim2,7:15:end)-ref1)./ref1];    % 1st configuration
ratio2=[sim(jsim2,1),100*(sim(jsim2,11:15:end)-ref2)./ref2];   % 2nd  configurationn
ratio3=[sim(jsim2,1),100*(sim(jsim2,13:15:end)-ref3)./ref3];   % 1st configuraiton + SL correction 

m=[sim(jsim2,1),sim(jsim2,4:15:end)]; % airmas

x1=smoothdata(ratio1(:,2:end),'gaussian',15,'SamplePoints',ratio1(:,1));
x2=smoothdata(ratio2(:,2:end),'gaussian',15,'SamplePoints',ratio2(:,1));
x3=smoothdata(ratio3(:,2:end),'gaussian',15,'SamplePoints',ratio3(:,1));

f=figure;
set(f,'Tag','comp_2019');
plot(ratio1(:,1),x1,'.'); hold all
grid;
datetick('x',12);
title(strcat('Conf. ',confname(j)))
set(gca,'Ylim',[-2,2]);
legend(string(Cal.brw));
vline_v(bevents.Date_mat,'-',bevents.Row)



%tabla por osc
m=[sim(jsim2,1),sim(jsim2,4:15:end)];
osc=[m(:,1),ref1.*m(:,2:end)];
r_osc=[ratio1,osc(:,2:end)];

bevents=brw_events(brw_events.Date>=datetime(2019,05,24),:)
[g,p]=group_time_new(r_osc,bevents.Date_mat);
gr=unique(g);
if any(gr==0)
   gr(gr==0)=[];
   Nobs=sum(p(:,gr+1))
   gr=gr+1;
end

p=logical(p);

for jj=1:length(gr)
   [r1,b,c]=group_var(r_osc(p(:,gr(jj)),:),[200,400,700,1000,1500]);
   [m,s]=grpstats(r1,r1(:,end));

t_osc{jj,j}=array2table(m(:,2:5),'RowNames',c(m(:,end)),'VariableNames',varname([Cal.brw_name,{'osc'}]));
t_osc{jj,j}.Properties.Description = bevents.Row{gr(jj)-1};

end

end


%% Comparamos solo 157 y 183
for j=1:2
 
TSYNC=10; % (min)
ref1=fix(ds{1,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
ref2=fix(ds{2,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
s=scan_join([ref1,ds{1,j}],[ref2,ds{2,j}]);
label{1}='Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc';
% s 1 date
%        2:16  157 'Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc'
%       17:31  183
%       32:45  185

% datos simultaneos
jsim=all(~isnan(s(:,7:15:end)),2) ; 
ref=median(s(jsim,7:15:end),2);
sim=s(jsim,:);

jsim2=all(abs(sim(:,3:15:end)-mean(sim(:,3:15:end),2))<1.0,2);

ref1=median(sim(jsim2,7:15:end),2);  % reference 1st configuration  (O3_1)
ref2=median(sim(jsim2,11:15:end),2);  % reference 2st configuration  (O3_2)
ref3=median(sim(jsim2,13:15:end),2);  % reference 1st configuration+SL  (O3_1_sl)

ratio1=[sim(jsim2,1),100*(sim(jsim2,7:15:end)-ref1)./ref1];    % 1st configuration
ratio2=[sim(jsim2,1),100*(sim(jsim2,11:15:end)-ref2)./ref2];   % 2nd  configurationn
ratio3=[sim(jsim2,1),100*(sim(jsim2,13:15:end)-ref3)./ref3];   % 1st configuraiton + SL correction 

m=[sim(jsim2,1),sim(jsim2,4:15:end)]; % airmas
x1=smoothdata(ratio1(:,2:end),'gaussian',15,'SamplePoints',ratio1(:,1));
x2=smoothdata(ratio2(:,2:end),'gaussian',15,'SamplePoints',ratio2(:,1));
x3=smoothdata(ratio3(:,2:end),'gaussian',15,'SamplePoints',ratio3(:,1));

f=figure;
set(f,'Tag','comp_2019');
plot(ratio1(:,1),x1,'.'); hold all
%plot(ratio2(:,1),x2,'-','LineWidth',2); hold all
%plot(ratio2(:,1),x3,':'); hold all
grid;
datetick('x',12);
%title(['Configuracion',num2str(j)])
title(strcat('Conf. ',confname(j)))
set(gca,'Ylim',[-2,2]);
legend(string(Cal.brw));
% for ei=1:2
% vline_v(events{ei}.dates,'k',strrep(events{ei}.labels,'_','\_'));
% end
vline_v(bevents.Date_mat,'-',bevents.Row)
end


ix=sort(findobj('-regexp','Tag','comp_2019'));
Width=20; Height=10;
printfiles_report(ix',Cal.dir_figs,'Width',Width,'Height',Height);

%tabla por osc
m=[sim(jsim2,1),sim(jsim2,4:15:end)];
osc=[m(:,1),ref1.*m(:,2:end)];
r_osc=[ratio1,osc(:,2:end)];

bevents=brw_events(brw_events.Date>=datetime(2019,05,24),:)
[g,p]=group_time(r_osc,bevents.Date_mat);
[r1,b,c]=group_var(r_osc,[200,400,700,1000,1500])
gr=find(unique(g));
p=logical(p);
%% dataset brw, time groups, data
dataset=NaN*ones(2,length(gr),length(r_osc));
for bb=1:2
 for jj=1:length(gr)
    dataset(bb,jj,p(:,gr(jj)))=ratio1(p(:,gr(jj)),bb+1);
 end
end
figure
h=boxplot2(dataset,'whisker',3)
cmap = get(0, 'defaultaxescolororder');
for ii = 1:3
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h);
end
set([h.lwhis h.uwhis], 'linestyle', '-','LineWidth',2);
set(h.out, 'marker', '.');
set(gca,'XtickLabels',[{''},labels',{''},{''},labels',{''}])
xtickangle(30)
title({'Brewer #157 (left) and Brewer#183 (right) ','Before , during and affer Huelva campaing'})
ylabel(' % difference ratio')
grid on

%% dataset
%% dataset brw, time groups, data
%dataset=NaN*ones(2,length(gr),5,length(r_osc));
osc_range=[400,500,700,1000];
dataset=NaN*ones(2,length(gr),length(osc_range),length(r_osc));

%
for jj=1:length(gr)
   [r1,b,c]=group_var(r_osc(p(:,gr(jj)),:),osc_range);
   [m,s]=grpstats(r1,r1(:,end));
   for k=1:length(c)
    dataset(1,jj,k,r1(:,end)==k)=r1(r1(:,end)==k,2);
    dataset(2,jj,k,r1(:,end)==k)=r1(r1(:,end)==k,3);
   end

t1_osc{jj}=array2table(m(:,2:4),'RowNames',c(m(:,end)),'VariableNames',varname({'B#157','B#183','osc'}));
t1_osc{jj}.Properties.Description = bevents.Row{jj};

end
%
figure;h(1)=boxplot2(squeeze(dataset(1,:,:,:)),'whisker',3,'notch','on');
hold on;h(2)=boxplot2(squeeze(dataset(2,:,:,:)),'whisker',3,'notch','on');
for jj=1:2 % brewer
 for ii = 1:size(h(1).box,1)
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h(jj));
   
 end
%set([h(jj).lwhis h(jj).uwhis], 'linestyle', '-','LineWidth',2*jj);
set([h(jj).med], 'linestyle', '-','LineWidth',2*jj);
set(h(jj).out, 'marker', '.');
end
%set(gca,'XtickLabels',[{''},labels',{''},{''},labels
set(gca,'XtickLabels',[labels])
set(gca,'Xtick',1:3)
%xtickangle(30)
title({'Brewer #157 (thick line) and Brewer#183 ratio at different osc ranges ','Before , during and affer Huelva campaing'})
ylabel(' % difference ratio')
grid on
legend(h(1).med(:,1),c,'orientation','horizontal')

%%